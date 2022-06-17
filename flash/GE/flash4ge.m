function flash4ge(sys, N, FOV, FLIP, TR, DTE)
%
% Fully-sampled 3D RF-spoiled GRE sequence.
% Combines several acquisitions with different flip/tr/te.
% 
% Inputs:
%  sys        struct   scanner hardware settings. See toppe.systemspecs().
%  N          [1 3]    matrix size
%  FOV        [1 3]    field of view (cm)
%  FLIP       [1 n]    flip angle (degrees) for each of n scans
%  TR         [1 n]    sequence repetition time (ms) for the n scans
%  DTE        [1 n]    increase (from minimum value ) in TE for each of the n scans

nScans = numel(DTE);

if numel(FLIP) ~= nScans | numel(TR) ~= nScans
    error('FLIP/TR/DTE must be the same length');
end

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
ex.dur = 2;           % pulse duration (ms)
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

% Write entry file.
% This can be edited by hand as needed after copying to scanner.
fid = fopen('toppeN.entry', 'wt');
fprintf(fid, '/usr/g/research/pulseq/cal/\n');  
fprintf(fid, 'modules.txt\n');
fprintf(fid, 'scanloop.txt\n');
fprintf(fid, '%s\n', mods.ex);
fprintf(fid, '%s\n', mods.readout);
fprintf(fid, 'seqstamp.txt');
fclose(fid);


%% Create .mod files

% excitation module
[ex.rf, ex.g] = toppe.utils.rf.makeslr(max(FLIP), ex.slThick, ...
    ex.tbw, ex.dur, nz*nCyclesSpoil, sys, ...
    'ftype', ex.ftype, ...
    'spoilDerate', 0.5, ...
    'ofname', mods.ex);

% readout module
% Here we use the helper function 'makegre' to do that, but
% that's not a requirement.
% Reduce slew to keep PNS in normal mode (<80% of limit)
toppe.utils.makegre(FOV(1), N(1), voxSize(3), sys, ... 
    'ofname', mods.readout, ...
    'ncycles', nCyclesSpoil); 

% Display .mod files.
%toppe.plotmod(mods.ex, 'gradcoil', sys.gradient);
%toppe.plotmod(mods.readout, 'gradcoil', sys.gradient);


%% Get minimum TR
% We do this by building one TR
toppe.write2loop('setup', sys);
toppe.write2loop(mods.ex, sys);
toppe.write2loop(mods.readout, sys);
toppe.write2loop('finish', sys);

trmin = toppe.getTRtime(1, 2, sys)*1e3;  % ms


%% Write scanloop.txt
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;    % degrees
ny = N(2);
nz = N(3);

toppe.write2loop('setup', sys, 'version', 4);  % initialize file ('scanloop.txt')

for iim = 1:nScans  % loop over the various 3D image acquisitions
    for ii = 1:30
         fprintf('\b');
    end
    fprintf('Scan %d of %d', iim, nScans);

    rfamp = FLIP(iim)/max(FLIP);        % RF amplitude scaling
    dte = DTE(iim);                     % insert delay after RF pulse (delay echo)
    textra = max(0, TR(iim) - dte - trmin);   % delay after readout module

    % partition-encoding loop
    for iz = -2:nz     % We use iz<1 for approach to steady-state

        % phase-encoding loop
        for iy = 1:ny
            % y/z phase encode amplitudes. Turn off during approach to steady-state.
            % Convention: start at (-kymax, -kzmax)
            a_gy = -((iy-1+0.5)-ny/2)/(ny/2) * (iz>0);  
            a_gz = -((iz-1+0.5)-nz/2)/(nz/2) * (iz>0);

            toppe.write2loop(mods.ex, sys, ...
                'RFamplitude', rfamp, ...
                'textra', dte, ...
                'RFphase', rfphs);

            toppe.write2loop(mods.readout, sys, ...
                'Gamplitude', [1.0 a_gy a_gz]', ...
                'DAQphase', rfphs, ...
                'textra', textra, ...
                'slice', max(iz,1), 'echo', iim, 'view', iy);

            % update rf/rec phase
            rfphs = rfphs + (rf_spoil_seed/180*pi)*rf_spoil_seed_cnt ;  % radians
            rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
        end
    end
end
fprintf('\n');
toppe.write2loop('finish', sys);  % finalize file


%% Create 'sequence stamp' file for TOPPE
% This file is listed in line 6 of toppeN.entry
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', sys);


%% Write files to tar archive (for convenience).
system('tar cf flash.tar toppeN.entry seqstamp.txt scanloop.txt modules.txt *.mod');

toppe.utils.scanmsg('toppeN.entry');

% Play sequence in loop (movie) mode
%nModulesPerTR = 2;
%toppe.playseq(nModulesPerTR, sys, ...
%    'tpause', 0.05, ...
%    'nTRskip', 8);

return;

