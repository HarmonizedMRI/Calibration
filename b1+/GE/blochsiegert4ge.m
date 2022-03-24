function blochsiegert4ge(sys, N, FOV, nShots, tr, varargin)
% function blochsiegert4ge(sys, N, FOV, tr, [etl, Ry, isFlyback, isRampsamp, bsFreq)
%
% 3D FLASH segmented EPI (or spin warp) Bloch-Siegert B1+ mapping sequence.
%
% The following parameters are hard-coded:
%  Fermi function pulse shape of amplitude 0.05 Gauss
%  imaging flip angle 15 deg
%
% Inputs:
%  sys          struct   scanner hardware settings. See toppe.systemspecs().
%  N            [1 3]    matrix size
%  FOV          [1 3]    field of view (cm)
%  nShots       [1 1]    number of EPI shots (segments)
%  tr           [1 1]    sequence TR (ms)
% 
% Input options (keyword-value arguments):
%  bsFreq       [1 (>=2)]    Bloch-Siegert pulse frequency offsets (Hz). 
%                        Default: [-4000 4000]
%  Ry           [1 1]    EPI undersampling factor. Default: 1
%  flyback      bool     Flyback EPI? Default: true
%  rampsamp     bool     Ramp sampling? Default: false
%  decimation   int      Design EPI readout as if ADC dwell time is 4us*decimation.
%                        Default: 1
%                        Applies only to readouts w/o ramp sampling.
%                        (The actual ADC dwell time is fixed to 4us in TOPPE)

%% Parse input options
arg.bsFreq = [-4000 4000];  % Hz
arg.Ry = 1;
arg.flyback = true;
arg.rampsamp = false;
arg.decimation = 1;

% Substitute varargin values as appropriate and check inputs
arg = toppe.utils.vararg_pair(arg, varargin);


%% Set sequence parameters and perform other misc tasks

flip = 15;             % low flip angle for spin-density weighting
slThick = 0.8*FOV(3);  % slab thickness

% excitation (imaging) pulse
ex.tbw = 8;           % time-bandwidth product of SLR pulse 
ex.dur = 3;           % pulse duration (ms)
ex.ftype = 'min';     

% Fermi pulse
bs.amp = 0.05;            % Amplitude of Fermi pulse (Gauss)

nScans = numel(arg.bsFreq);

ny = N(2);
nz = N(3);

% .mod file names
mods.ex         = 'ex.mod';    % does not contain the rephasing gradient lobe
mods.bs         = 'bs.mod';    % BS pulse (inserted between ex.mod and ex-rephaser.mod)
mods.exRephaser = 'ex-rephaser.mod';  % slice-select rephaser
mods.prephaser  = 'prephaser.mod';  % for EPI readout
mods.readout    = 'readout.mod';    % EPI echo-train

fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', length(fieldnames(mods)));
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\n', mods.ex);
fprintf(fid, '%s\t0\t1\t0\n', mods.bs);
fprintf(fid, '%s\t0\t0\t0\n', mods.exRephaser);
fprintf(fid, '%s\t0\t0\t0\n', mods.prephaser);
fprintf(fid, '%s\t0\t0\t1\n', mods.readout);
fclose(fid);


%% Create .mod files 

% readout.mod and prephaser.mod
[gx,gy,gz] = toppe.utils.makeepi(FOV, N, ...
    nShots, sys, ...
    'flyback', arg.flyback, ...
    'rampsamp', arg.rampsamp, ...
    'Ry', arg.Ry, ...
    'writefiles', true); 

% echo spacing (ms)
tsamp = sys.raster*1e3;   % gradient sample time (ms)
if arg.flyback
    es = tsamp * (numel(gx.echo) + numel(gx.flyback));
else
    es = tsamp * numel(gx.echo);
end
if nShots == N(2)
    es = 0;
end

% imaging excitation pulse
nCyclesSpoil = 4*nz;  % unbalanced gradients
[ex.rf, ex.g] = toppe.utils.rf.makeslr(flip, slThick, ...
        ex.tbw, ex.dur, nCyclesSpoil, sys, ...
        'ftype', ex.ftype, ...
        'spoilDerate', 0.5, ...
        'writeModFile', false);

% Split off gradient rephaser lobe, and write to separate .mod files
I = 4 + find(ex.g(5:end) < 0);
irep = I(1);    % start of rephaser lobe
ex.gex = [ex.g(1:(irep-1))];
ex.rephaser = [0; ex.g(irep:end)];

ex.rf = toppe.makeGElength(ex.rf(1:length(ex.gex)));
ex.gex = toppe.makeGElength(ex.gex);
ex.rephaser = toppe.makeGElength(ex.rephaser);

toppe.writemod(sys, 'rf', ex.rf, 'gz', ex.gex, ...
    'ofname', mods.ex);

toppe.writemod(sys, 'gz', ex.rephaser, ...
    'ofname', mods.exRephaser);

% Bloch-Siegert module
toppe.utils.rf.makebs(bs.amp, 'system', sys, 'ofname', mods.bs);


%% Get minimum TR
% We do this by building one TR
toppe.write2loop('setup', sys);
toppe.write2loop(mods.ex, sys);
toppe.write2loop(mods.bs, sys);
toppe.write2loop(mods.exRephaser, sys);
toppe.write2loop(mods.prephaser, sys);
toppe.write2loop(mods.readout, sys);
toppe.write2loop(mods.prephaser, sys);
toppe.write2loop('finish', sys);

trmin = (nShots-1)/nShots*es + toppe.getTRtime(1, 6, sys)*1e3;  % ms


%% Write scanloop.txt

rfphs = 0;    % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

toppe.write2loop('setup', sys);
for iim = 1:nScans
    for iz = 0:nz   % iz < 1 are discarded acquisitions to reach steady state
        for iy = 1:nShots
            % y/z phase-encode amplitudes, scaled to (-1,1)
            a_gy = -(iz>0)*((iy-1+0.5)-ny/2)/(ny/2);
            a_gz = -(iz>0)*((iz-1+0.5)-nz/2)/(nz/2);   

            % rf excitation, bloch-siegert pulse, and slice rephaser
            toppe.write2loop(mods.ex, sys, ...
                'RFphase', rfphs);
            toppe.write2loop(mods.bs, sys, ...
                'RFoffset', arg.bsFreq(iim));
            toppe.write2loop(mods.exRephaser, sys, ...
                'textra', (iy-1)/nShots*es);

            % prephase (move to corner of kspace)
            toppe.write2loop(mods.prephaser, sys, ...
                'Gamplitude', [1 a_gy a_gz]');

            % readout. Data is stored in 'slice', 'echo', and 'view' indeces.
            toppe.write2loop(mods.readout, sys, ...
                'DAQphase', rfphs, ...
                'slice', iim, 'echo', max(iz,1), 'view', iy, ...
                'textra', es - iy/nShots*es, ...
                'dabmode', 'on');

            % rephase y and z encoding gradients (required condition for steady state)
            toppe.write2loop(mods.prephaser, sys, ...
                'textra', max(0, tr - trmin), ...
                'Gamplitude', [1 -a_gy -a_gz]');

            % update rf phase (RF spoiling)
            rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
            rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
        end
    end
end
fprintf('\n');
toppe.write2loop('finish', sys);

return

%% Create 'sequence stamp' file for TOPPE
% This file is listed in the 5th row in toppeN.entry
% NB! The file toppeN.entry must exist in the folder from where this script is called.
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', sys);


%% create tar file
system('tar czf ~/tmp/bs.tgz toppeN.entry seqstamp.txt scanloop.txt modules.txt *.mod');

return;

%% simulate slice profile
[rf,~,~,gz] = toppe.readmod('ex.mod');
X = linspace(-sp.rf.slThick, sp.rf.slThick, 100);         % cm
dt = 4e-3;            % ms
Gamp = max(gz);       % Gauss
T1 = 500;
T2 = 50;
toppe.utils.rf.slicesim(rf,dt,T1,T2,Gamp,X);

