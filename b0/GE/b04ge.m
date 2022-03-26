function b04ge(sys, N, FOV, flip, DTE, varargin)
% function b04ge(sys, N, FOV, flip, DTE, varargin)
%
% Fully-sampled 3D RF-spoiled GRE sequence for B0 (and B1-) mapping.
% 
% Inputs:
%  sys        struct   scanner hardware settings. See toppe.systemspecs().
%  N          [1 3]    matrix size
%  FOV        [1 3]    field of view (cm)
%  flip       [1 1]    flip angle (degrees)
%  DTE        [1 n]    increase (from minimum value ) in TE for each of n scans
%
% Input options with defaults:
%  entryFile = 'toppeN.entry';
%  filePath = '/usr/g/research/pulseq/cal/b0/';  % put scan files here
%  tbw = 8;                     % RF pulse time-bandwidth product
%  rfDur = 2;                   % RF pulse duration (ms)
%  ftype = 'min';               % 'min': minimum-phase SLR pulse; 'ls': linear phase
%  slabThick = 0.8*FOV(3);      % excited slab thickness
%  rfSpoilSeed = 117;         % RF spoiling phase increment factor (degrees)
%  exMod         = 'tipdown.mod';
%  readoutMod    = 'readout.mod';
%  nCyclesSpoil = 2;   % number of cycles of phase across voxel (along x and z)

% defaults
arg.entryFile = 'toppeN.entry';
arg.filePath = '/usr/g/research/pulseq/cal/b0/';
arg.tbw = 8;                     % RF pulse time-bandwidth product
arg.rfDur = 2;                   % RF pulse duration (ms)
arg.ftype = 'min';               % 'min': minimum-phase SLR pulse; 'ls': linear phase
arg.slabThick = 0.8*FOV(3);      % excited slab thickness
arg.rfSpoilSeed = 117;         % RF spoiling phase increment factor (degrees)
arg.exMod         = 'tipdown.mod';
arg.readoutMod    = 'readout.mod';
arg.nCyclesSpoil = 2;   % number of cycles of phase across voxel (along x and z)

% substitute with provided keyword arguments
arg = toppe.utils.vararg_pair(arg, varargin);

nScans = numel(DTE);

% Since we are using the helper function 'makegre' below,
% the in-plane FOV and matrix size must be square.
if N(1) ~= N(2) | FOV(1) ~= FOV(2)
    error('In-plane FOV and matrix be square.');
end

voxSize = FOV./N;  % cm

ny = N(2);
nz = N(3);

arg.exMod         = 'tipdown.mod';
arg.redaoutmod    = 'readout.mod';

% Write modules.txt
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', 2);
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\n', arg.exMod);
fprintf(fid, '%s\t0\t0\t1\n', arg.redaoutmod);
fclose(fid);

% Write entry file.
% This can be edited by hand as needed after copying to scanner.
toppe.writeentryfile(arg.entryFile, ...
    'filePath', arg.filePath, ...
    'b1ScalingFile', arg.exMod, ...
    'readoutFile', arg.redaoutmod);

%% Create .mod files

% excitation module
[ex.rf, ex.g] = toppe.utils.rf.makeslr(flip, arg.slabThick, ...
    arg.tbw, arg.rfDur, nz*arg.nCyclesSpoil, sys, ...
    'ftype', arg.ftype, ...
    'spoilDerate', 0.5, ...
    'ofname', arg.exMod);

% readout module
% Here we use the helper function 'makegre' to do that, but
% that's not a requirement.
% Reduce slew to keep PNS in normal mode (<80% of limit)
toppe.utils.makegre(FOV(1), N(1), voxSize(3), sys, ... 
    'ofname', arg.redaoutmod, ...
    'ncycles', arg.nCyclesSpoil); 


%% Write scanloop.txt
rfphs = 0;              % radians
rfSpoilSeed_cnt = 0;
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

            toppe.write2loop(arg.exMod, sys, ...
                'RFamplitude', 1.0, ...
                'textra', DTE(ite), ...
                'RFphase', rfphs);

            toppe.write2loop(arg.redaoutmod, sys, ...
                'Gamplitude', [1.0 a_gy a_gz]', ...
                'DAQphase', rfphs, ...
                'textra', max(DTE) - DTE(ite), ... % to keep TR constant
                'slice', max(iz,1), 'echo', ite, 'view', iy);

            % update rf/rec phase
            rfphs = rfphs + (arg.rfSpoilSeed/180*pi)*rfSpoilSeed_cnt ;  % radians
            rfSpoilSeed_cnt = rfSpoilSeed_cnt + 1;
        end
    end
end
fprintf('\n');
toppe.write2loop('finish', sys);  % finalize file

fprintf('TR = %.3f ms\n', toppe.getTRtime(1, 2, sys)*1e3);

%figure; toppe.plotseq(1, 4, sys);

% Create 'sequence stamp' file for TOPPE
% This file is listed in line 6 of the .entry file
toppe.preflightcheck(arg.entryFile, 'seqstamp.txt', sys);

% Write files to tar archive (for convenience).
system(sprintf('tar cf b0.tar %s seqstamp.txt scanloop.txt modules.txt *.mod', arg.entryFile));

toppe.utils.scanmsg(arg.entryFile);

% Play sequence in loop (movie) mode
%nModulesPerTR = 2;
%toppe.playseq(nModulesPerTR, sys, ...
%    'tpause', 0.05, ...
%    'nTRskip', 8);

return;

