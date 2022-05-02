function blochsiegert4ge(sys, N, FOV, nShots, tr, bsFreq, varargin)
% function blochsiegert4ge(sys, N, FOV, tr, [etl, Ry, isFlyback, isRampsamp, bsFreq)
%
% 3D FLASH segmented EPI (or spin warp) Bloch-Siegert B1+ mapping sequence.
%
% The following parameters are hard-coded (edit as needed):
%  Fermi function pulse shape of amplitude 0.05 Gauss
%  imaging flip angle 15 deg
%  .entry file name
%
% Inputs:
%  sys          struct   scanner hardware settings. See toppe.systemspecs().
%  N            [1 3]    matrix size
%  FOV          [1 3]    field of view (cm)
%  nShots       [1 1]    number of EPI shots (segments)
%  tr           [1 1]    sequence TR (ms)
%  bsFreq       [1 1]    Bloch-Siegert pulse frequency offset (Hz)
% 
% Input options (keyword-value arguments) with defaults:
%  entryFile  = 'toppeN.entry';
%  filePath   = '/usr/g/research/pulseq/cal/b1/'   % location of toppe<N>.entry file
%  Ry         = 1;         % EPI undersampling factor
%  flyback    = true;      % flyback EPI or not
%  rampsamp   = false;     % sample on ramps or not
%  decimation = 1;         % Design EPI readout as if ADC dwell time is 4us*decimation
                           % Applies only to readouts w/o ramp sampling.
                           % (The actual ADC dwell time is fixed to 4us in TOPPE)
%  flip = 15;                % degrees
%  slabThick = 0.8*FOV(3);   % slab thickness
%  tbw = 8;           % time-bandwidth product of SLR pulse 
%  rfdur = 3;         % RF pulse duration (ms)
%  ftype = 'min';     % minimum-phase (min) or linear phase (ls)
%  bsamp = 0.05;      % Amplitude of Fermi pulse (Gauss)
%
% Example usage:
%   sys = toppe.systemspecs('maxSlew', 10, 'gradient', 'xrm');                
%   N = [100 100 8];    % matrix size
%   FOV = [20 20 1.6];  % cm
%   nShots = 5;         % number of EPI segments
%   tr = 150;           % sequence TR (ms)
%   blochsiegert4ge(sys, N, FOV, nShots, tr, 'flyback', false, 'rampsamp', true);
%   toppe.plotseq(1,6,sys);


%% Parse input options
arg.entryFile = 'toppeN.entry';
arg.filePath = '/usr/g/research/pulseq/cal/b1/';
arg.Ry = 1;                 % EPI undersampling factor
arg.flyback = true;         % flyback EPI or not
arg.rampsamp = false;       % sample on ramps or not
arg.decimation = 1;         % Design EPI readout as if ADC dwell time is 4us*decimation
                            % Applies only to readouts w/o ramp sampling.
                            % (The actual ADC dwell time is fixed to 4us in TOPPE)
arg.flip = 15;                % degrees
arg.slabThick = 0.8*FOV(3);   % slab thickness
arg.tbw = 8;           % time-bandwidth product of SLR pulse 
arg.rfdur = 3;         % RF pulse duration (ms)
arg.ftype = 'min';     % minimum-phase (min) or linear phase (ls)
arg.bsamp = 0.05;      % Amplitude of Fermi pulse (Gauss)

% Substitute varargin values as appropriate and check inputs
arg = toppe.utils.vararg_pair(arg, varargin);


%% Set sequence parameters and perform other misc tasks

ny = N(2);
nz = N(3);

% .mod file names
mods.ex         = 'tipdown.mod';    % does not contain the rephasing gradient lobe
mods.bs         = 'bs.mod';    % BS pulse (inserted between tipdown.mod and tipdown-rephaser.mod)
mods.exRephaser = 'tipdown-rephaser.mod';  % slice-select rephaser
mods.prephaser  = 'prephaser.mod';  % for EPI readout
mods.readout    = 'readout.mod';    % EPI echo-train

% Write modules.txt
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

% Write entry file.
% This can be edited by hand as needed after copying to scanner.
toppe.writeentryfile(arg.entryFile, ...
    'filePath', arg.filePath, ...
    'b1ScalingFile', mods.ex, ...
    'readoutFile', mods.readout);


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
nCyclesSpoil = 2*nz;  % unbalanced gradients
[ex.rf, ex.g] = toppe.utils.rf.makeslr(arg.flip, arg.slabThick, ...
        arg.tbw, arg.rfdur, nCyclesSpoil, sys, ...
        'ftype', arg.ftype, ...
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
toppe.utils.rf.makebs(arg.bsamp, 'system', sys, 'ofname', mods.bs);


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

freq = [-bsFreq bsFreq];

toppe.write2loop('setup', sys);
for iz = 0:nz   % iz < 1 are discarded acquisitions to reach steady state
    for iy = 1:nShots
        for iim = 1:length(freq)
            % y/z phase-encode amplitudes, scaled to (-1,1)
            a_gy = -(iz>0)*((iy-1+0.5)-ny/2)/(ny/2);
            a_gz = -(iz>0)*((iz-1+0.5)-nz/2)/(nz/2);   

            % rf excitation, bloch-siegert pulse, and slice rephaser
            toppe.write2loop(mods.ex, sys, ...
                'RFphase', rfphs);
            toppe.write2loop(mods.bs, sys, ...
                'RFoffset', freq(iim));
            toppe.write2loop(mods.exRephaser, sys, ...
                'textra', (iy-1)/nShots*es);

            % prephase (move to corner of kspace)
            toppe.write2loop(mods.prephaser, sys, ...
                'Gamplitude', [1 a_gy a_gz]');

            % readout. Data is stored in 'slice', 'echo', and 'view' indeces.
            toppe.write2loop(mods.readout, sys, ...
                'DAQphase', rfphs, ...
                'echo', iim, 'slice', max(iz,1), 'view', iy, ...
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

% Create 'sequence stamp' file for TOPPE
% This file is listed in the 5th row in toppeN.entry
toppe.preflightcheck(arg.entryFile, 'seqstamp.txt', sys);

% create tar file (for convenience)
system(sprintf('tar cf b1.tar %s seqstamp.txt scanloop.txt modules.txt *.mod', arg.entryFile));

toppe.utils.scanmsg(arg.entryFile);

return;

%% simulate slice profile
[rf,~,~,gz] = toppe.readmod('ex.mod');
X = linspace(-sp.rf.slThick, sp.rf.slThick, 100);         % cm
dt = 4e-3;            % ms
Gamp = max(gz);       % Gauss
T1 = 500;
T2 = 50;
toppe.utils.rf.slicesim(rf,dt,T1,T2,Gamp,X);

