function b04ge(sys, N, FOV, flip, DTE)
% function b04ge(sys, N, FOV, flip, DTE)
%
% Fully-sampled 3D RF-spoiled GRE sequence for B0 (and B1-) mapping.
% 
% Inputs:
%  sys        struct   scanner hardware settings. See toppe.systemspecs().
%  N          [1 3]    matrix size
%  FOV        [1 3]    field of view (cm)
%  flip       [1 1]    flip angle (degrees)
%  DTE        [1 n]    increase (from minimum value ) in TE for each of n scans

nScans = numel(DTE);

% Since we are using the helper function 'makegre' below,
% the in-plane FOV and matrix size must be square.
if N(1) ~= N(2) | FOV(1) ~= FOV(2)
    error('In-plane FOV and matrix be square.');
end

% number of cycles of spoiler phase across voxel dimension (applied along x and z)
nCyclesSpoil = 2;           

voxSize = FOV./N;  % cm

ny = N(2);
nz = N(3);

% excitation (imaging) pulse parameters
ex.tbw = 8;           % time-bandwidth product of SLR pulse 
ex.dur = 3;           % pulse duration (ms)
ex.ftype = 'min';     % minimum-phase SLR pulse
ex.slThick = 0.8*FOV(3);  % slab thickness

% .mod file names
mods.ex         = 'tipdown.mod';
mods.readout    = 'readout.mod';

% Write modules.txt
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', 2);
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\n', mods.ex);
fprintf(fid, '%s\t0\t0\t1\n', mods.readout);
fclose(fid);


%% Create .mod files

% excitation module
[ex.rf, ex.g] = toppe.utils.rf.makeslr(flip, ex.slThick, ...
        ex.tbw, ex.dur, nz*nCyclesSpoil, sys, ...
        'ftype', ex.ftype, ...
        'spoilDerate', 0.5, ...
        'ofname', mods.ex);

% readout module
% Here we use the helper function 'makegre' to do that, but
% that's not a requirement.
% Reduce slew to keep PNS in normal mode (<80% of limit)
systmp = sys;
systmp.maxSlew = 15;
toppe.utils.makegre(FOV(1), N(1), voxSize(3), systmp, ... 
    'ofname', mods.readout, ...
    'ncycles', nCyclesSpoil); 

% Display .mod files.
toppe.plotmod(mods.ex, 'gradcoil', sys.gradient);
toppe.plotmod(mods.readout, 'gradcoil', sys.gradient);


%% Write scanloop.txt
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;    % degrees
ny = N(2);
nz = N(3);

toppe.write2loop('setup', sys, 'version', 4);  % initialize file ('scanloop.txt')

for iz = -1:nz     % We use iz<1 for approach to steady-state
    fprintf('\b\b\b\b\b\b\b\b%d of %d', max(1,iz), nz);
    for iy = 1:ny
        for ite = 1:length(DTE)
            % y/z phase encode amplitudes. Turn off during approach to steady-state.
            % My convention is to start at (-kymax, -kzmax)
            a_gy = -((iy-1+0.5)-ny/2)/(ny/2) * (iz>0);  
            a_gz = -((iz-1+0.5)-nz/2)/(nz/2) * (iz>0);

            toppe.write2loop(mods.ex, sys, ...
                'RFamplitude', 1.0, ...
                'textra', DTE(ite), ...
                'RFphase', rfphs);

            toppe.write2loop(mods.readout, sys, ...
                'Gamplitude', [1.0 a_gy a_gz]', ...
                'DAQphase', rfphs, ...
                'textra', max(DTE) - DTE(ite), ... % to keep TR constant
                'slice', max(iz,1), 'echo', ite, 'view', iy);

            % update rf/rec phase
            rfphs = rfphs + (rf_spoil_seed/180*pi)*rf_spoil_seed_cnt ;  % radians
            rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
        end
    end
end
fprintf('\n');
toppe.write2loop('finish', sys);  % finalize file

fprintf('TR = %.3f ms\n', toppe.getTRtime(1, 2, sys)*1e3);

figure; toppe.plotseq(1, 4, sys);

return


%% Create 'sequence stamp' file for TOPPE
% This file is listed in line 6 of toppeN.entry
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', sys);


%% Write files to tar archive (for convenience).
system('tar cf b0.tar toppeN.entry modules.txt scanloop.txt *.mod seqstamp.txt');

% Play sequence in loop (movie) mode
%nModulesPerTR = 2;
%toppe.playseq(nModulesPerTR, sys, ...
%    'tpause', 0.05, ...
%    'nTRskip', 8);

return;

