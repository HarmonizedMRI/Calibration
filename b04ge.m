function b04ge(sys)
% function b04ge(sys)
%
% Create fully-sampled 3D RF-spoiled GRE sequence for B1- and B0 mapping.
% 
% Usage:
%  1. Set hardware/design limits in getsys.m
%  2. >> sys = getsys;
%  3. >> b04ge(sys.ge);
%
% This script generates the following files:
%  modules.txt
%  scanloop.txt
%  tipdown.mod
%  readout.mod
%  seqstamp.txt

%% Acquisition parameters
imSize = [64 64 64];        % matrix size
FOV = [24 24 24];           % cm
flip = 3;                   % excitation flip angle (degrees). Low flip for m0 weighting
nCyclesSpoil = 2;           % number of cycles of spoiler phase across voxel dimension (applied along x and z)
deltaTE = 1000/440*[0 1];   % multiple TE times for B0 mapping (ms)

voxSize = FOV./imSize;  % cm

% Since we are using the helper function 'makegre' below,
% the in-plane FOV and matrix size must be square.
if imSize(1) ~= imSize(2) | FOV(1) ~= FOV(2)
    error('In-plane FOV and matrix be square.');
end


%% Create non-selective excitation module
% Use bandwidth-limited RF shape to mitigate poor RF system performance (??)
% Hard pulse with rounded edges (Fermi function, like in Bloch-Siegert mapping).
% Include spoiler gradient along z (more time efficient).
dur = 0.1;  % ms
width = 0.001e-3;  % Fermi transition width (s)
f = toppe.utils.rf.fermi(1, 0, dur/2*1e-3, width);
rf = (flip/360) / (sys.ge.gamma * sys.ge.raster * sum(f)) * f(:);  % Gauss
nChop = [0 8];  % to allow time for RF ringdown when converting to Pulseq
rf = [rf; zeros(nChop(2),1)]; 

if false
% Do Bloch simulations to check flip angle and frequency response
dt = sys.ge.raster*1e3;  % ms
m0 = [0 0 1];
Z = linspace(-5,5,100);
gz = 0*rf;
figure; toppe.utils.rf.slicesim(m0,rf,0*rf,dt,Z,1000,100);
figure; toppe.utils.rf.plotrfspectrum(rf, 'bo');
end

tmp = toppe.utils.makecrusher(nCyclesSpoil, voxSize(3), sys.ge, 0, sys.ge.maxSlew/2, sys.ge.maxGrad);
tmp = [tmp; zeros(5,1)];  % to make sure gradients are off during RF transmission
gSpoil = toppe.makeGElength([tmp; 0*rf]);
rf = toppe.makeGElength([0*tmp; rf]);
toppe.writemod(sys.ge, ...
    'rf', rf, 'gz', gSpoil, ...
    'nomflip', flip, ...
    'nChop', nChop, ...
    'ofname', 'tipdown.mod');


%% Create readout module
% Here we use the helper function 'makegre' to do that, but
% that's not a requirement.
% Reduce slew to keep PNS in normal mode (<80% of limit)
systmp = sys.ge;
systmp.maxSlew = 15;
toppe.utils.makegre(FOV(1), imSize(1), FOV(3)/imSize(3), systmp, ... 
    'ofname', 'readout.mod', ...
    'ncycles', nCyclesSpoil); 

% Display .mod files.
%toppe.plotmod('tipdown.mod', 'gradcoil', sys.ge.gradient);
%toppe.plotmod('readout.mod', 'gradcoil', sys.ge.gradient);


%% Write modules.txt
fid = fopen('modules.txt', 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', 2);
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\n', 'tipdown.mod');
fprintf(fid, '%s\t0\t0\t1\n', 'readout.mod');
fclose(fid);


%% Write scanloop.txt
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;    % degrees
ny = imSize(2);
nz = imSize(3);

toppe.write2loop('setup', sys.ge, 'version', 4);  % initialize file ('scanloop.txt')

for iz = -1:nz     % We use iz<1 for approach to steady-state
    fprintf('\b\b\b\b\b\b\b\b%d of %d', max(1,iz), nz);
    for iy = 1:ny
        for ite = 1:length(deltaTE)
            % y/z phase encode amplitudes. Turn off during approach to steady-state.
            % My convention is to start at (-kymax, -kzmax)
            a_gy = -((iy-1+0.5)-ny/2)/(ny/2) * (iz>0);  
            a_gz = -((iz-1+0.5)-nz/2)/(nz/2) * (iz>0);

            toppe.write2loop('tipdown.mod', sys.ge, ...
                'RFamplitude', 1.0, ...
                'textra', deltaTE(ite), ...
                'RFphase', rfphs);

            toppe.write2loop('readout.mod', sys.ge, ...
                'Gamplitude', [1.0 a_gy a_gz]', ...
                'DAQphase', rfphs, ...
                'textra', max(deltaTE) - deltaTE(ite), ... % to keep TR constant
                'slice', max(iz,1), 'echo', ite, 'view', iy);

            % update rf/rec phase
            rfphs = rfphs + (rf_spoil_seed/180*pi)*rf_spoil_seed_cnt ;  % radians
            rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
        end
    end
end
fprintf('\n');
toppe.write2loop('finish', sys.ge);  % finalize file

fprintf('TR = %.3f ms\n', toppe.getTRtime(1, 2, sys.ge)*1e3);


%% Create 'sequence stamp' file for TOPPE
% This file is listed in line 6 of toppeN.entry
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', sys.ge);


%% Write files to tar archive (for convenience).
system('tar cf b0.tar toppeN.entry modules.txt scanloop.txt *.mod seqstamp.txt');

% Play sequence in loop (movie) mode
%nModulesPerTR = 2;
%toppe.playseq(nModulesPerTR, sys.ge, ...
%    'tpause', 0.05, ...
%    'nTRskip', 8);

return;

