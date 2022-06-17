function [b1, rss] = loadblochsiegert4ge(pfile, readoutFile, bsModuleFile, bsFreq)
% function [b1, rss] = loadblochsiegert4ge(pfile, readoutFile)
%
% WIP. Not tested.
% TODO:
%  if non-flyback: do odd/even EPI echo correction
%  if ramp sampling: do interpolation
%  Add BS pulse info to bs.mod header, and load here
%
% Get b1+ map obtained with the TOPPE sequence files created
% with blochsiegert4ge.m
%
% Inputs:
%  pfile         string
%  readoutFile   string   Module file name, e.g., 'readout.mod'
%  bsMofuleFile  string   Module file name, e.g., 'bs.mod'
%  bsFreq        [1 1]    Bloch-Siegert pulse frequency offsets (Hz)
%
% Outputs:
%  b1        [nx ny nz]   b1 map (Gauss)
%  rss       [nx ny nz]   Root-sum-of-squares coil-combined image from first echo

% Get decimation from readout module file header
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod(readoutFile);
nx = paramsint16(1);   % image size
ny = paramsint16(2);
%Ry = paramsint16(4);
decimation = paramsint16(7);

% Load EPI raw data
for echo = 1:2
    d(:,:,:,:,echo) = toppe.utils.loadepi(pfile, echo, readoutFile);  % [nx*decimation ny nz nCoils]
end

% Reconstruct complex coil images
imNegFreq = toppe.utils.ift3(d(:,:,:,:,1));  % negative Bloch-Siegert frequency
imPosFreq = toppe.utils.ift3(d(:,:,:,:,2));  % 
bs.amp = max(abs(rf));           % Amplitude of Fermi pulse (Gauss)
rf = toppe.readmod(bsModuleFile);
dt = 4e-6;  % RF raster time (sec)
kbs = toppe.utils.rf.calckbs(rf, bsFreq, dt);
pc = toppe.utils.phasecontrast3d(imPosFreq, imNegFreq);
b1 = sqrt(pc/2/kbs);   % Measured b1 in Gauss. Note factor of 2!

