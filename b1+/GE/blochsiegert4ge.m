function b1(seq)
% Create segmented EPI (or spin warp) bSSFP scans

%% Get sequence parameters and define file names

nscans = numel(seq.bs.freq);

nx = seq.matrix(1);
ny = seq.matrix(2);
nz = seq.matrix(3);

moduleListFile = 'modules.txt';
% loopFile = 'scanloop.txt'; % this is the default name, so no need to define it
mods.ex         = 'ex.mod';    % does not contain the rephasing gradient lobe
mods.bs         = 'bs.mod';
mods.exRephaser = 'ex-rephaser.mod';
mods.readout    = seq.mods.readout;
mods.prephaser  = seq.mods.prephaser;


%% Write modules.txt
fid = fopen(moduleListFile, 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', length(fieldnames(mods)));
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\n', mods.ex);
fprintf(fid, '%s\t0\t1\t0\n', mods.bs);
fprintf(fid, '%s\t0\t0\t0\n', mods.exRephaser);
fprintf(fid, '%s\t0\t0\t1\n', mods.readout);
fprintf(fid, '%s\t0\t0\t0\n', mods.prephaser);
fclose(fid);


%% Create .mod files 
%% Create imaging and excitation modules (we'll reuse readout.mod from the bssfp scan)
% imaging excitation pulse
nCyclesSpoil = 4*nz;  % unbalanced gradients
[ex.rf, ex.g] = toppe.utils.rf.makeslr(seq.bs.flip, seq.rf.slThick, ...
        seq.rf.tbw, seq.rf.dur, nCyclesSpoil, seq.sys, ...
        'ftype', seq.rf.ftype, ...
        'spoilDerate', 0.5, ...
        'sliceOffset', seq.sliceOffset, ...
        'writeModFile', false);

% Split off gradient rephaser lobe,
% and write to separate .mod files
I = 4 + find(ex.g(5:end) < 0);
irep = I(1);    % start of rephaser lobe
ex.gex = [ex.g(1:(irep-1))];
ex.rephaser = [0; ex.g(irep:end)];

ex.rf = toppe.makeGElength(ex.rf(1:length(ex.gex)));
ex.gex = toppe.makeGElength(ex.gex);
ex.rephaser = toppe.makeGElength(ex.rephaser);

toppe.writemod(seq.sys, 'rf', ex.rf, 'gz', ex.gex, ...
    'ofname', mods.ex);

toppe.writemod(seq.sys, 'gz', ex.rephaser, ...
    'ofname', mods.exRephaser);

% Bloch-Siegert module
toppe.utils.rf.makebs(seq.bs.amp, 'system', seq.sys, 'ofname', mods.bs);


%% Get min TR so we can calculate textra (to achieve desired TR)
% We do this by building one TR
toppe.write2loop('setup', seq.sys);
toppe.write2loop(mods.ex, seq.sys);
toppe.write2loop(mods.bs, seq.sys);
toppe.write2loop(mods.exRephaser, seq.sys);
toppe.write2loop(mods.prephaser, seq.sys);
toppe.write2loop(mods.readout, seq.sys);
toppe.write2loop(mods.prephaser, seq.sys);
toppe.write2loop('finish', seq.sys);

trmin = (seq.nshots-1)/seq.nshots*seq.es + toppe.getTRtime(1, 6, seq.sys)*1e3;  % ms


%% Write scanloop.txt
ny = seq.matrix(2);
nz = seq.matrix(3);

rfphs = 0;    % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

toppe.write2loop('setup', seq.sys);
for iim = 1:nscans
    for iz = 0:nz   % iz < 1 are discarded acquisitions to reach steady state
        for iy = 1:seq.nshots
            % y/z phase-encode amplitudes, scaled to (-1,1)
            a_gy = -(iz>0)*((iy-1+0.5)-ny/2)/(ny/2);
            a_gz = -(iz>0)*((iz-1+0.5)-nz/2)/(nz/2);   

            % rf excitation, bloch-siegert pulse, and slice rephaser
            toppe.write2loop(mods.ex, seq.sys, ...
                'RFphase', rfphs);
            toppe.write2loop(mods.bs, seq.sys, ...
                'RFoffset', seq.bs.freq(iim));
            toppe.write2loop(mods.exRephaser, seq.sys, ...
                'textra', (iy-1)/seq.nshots*seq.es);

            % prephase (move to corner of kspace)
            toppe.write2loop(mods.prephaser, seq.sys, ...
                'Gamplitude', [1 a_gy a_gz]');

            % readout. Data is stored in 'slice', 'echo', and 'view' indeces.
            toppe.write2loop(mods.readout, seq.sys, ...
                'DAQphase', rfphs, ...
                'slice', iim, 'echo', max(iz,1), 'view', iy, ...
                'textra', seq.es - iy/seq.nshots*seq.es, ...
                'dabmode', 'on');

            % rephase y and z encoding gradients (required condition for steady state)
            toppe.write2loop(mods.prephaser, seq.sys, ...
                'textra', max(0, seq.bs.tr - trmin), ...
                'Gamplitude', [1 -a_gy -a_gz]');

            % update rf phase (RF spoiling)
            rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
            rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
        end
    end
end
fprintf('\n');
toppe.write2loop('finish', seq.sys);


%% Create 'sequence stamp' file for TOPPE.
% This file is listed in the 5th row in toppeN.entry
% NB! The file toppeN.entry must exist in the folder from where this script is called.
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', seq.sys);


%% create tar file
system('tar czf ~/tmp/scan,exchange,b1.tgz seqstamp.txt modules.txt scanloop.txt *.mod getparams.m b1.m');

return;

%% simulate slice profile
[rf,~,~,gz] = toppe.readmod('ex.mod');
X = linspace(-seq.rf.slThick, seq.rf.slThick, 100);         % cm
dt = 4e-3;            % ms
Gamp = max(gz);       % Gauss
T1 = 500;
T2 = 50;
toppe.utils.rf.slicesim(rf,dt,T1,T2,Gamp,X);

