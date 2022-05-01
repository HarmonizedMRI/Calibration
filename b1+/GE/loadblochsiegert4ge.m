function [b1, rss] = loadblochsiegert4ge(pfile, readoutFile, bsModuleFile)
% function [b1, rss] = loadblochsiegert4ge(pfile, readoutFile)
%
% Reconstruct individual coil images and B0 map
%
% Inputs:
%  pfile         string
%  readoutFile   string. Default: 'readout.mod' inside 'b0.tar'
%
% Outputs:
%  rssim         Root-sum-of-squares coil-combined image from first echo

[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod(readoutFile);
nx = paramsint16(1);   % image size
ny = paramsint16(2);
%Ry = paramsint16(4);
decimation = paramsint16(7);

for echo = 1:2
    d(:,:,:,:,echo) = toppe.utils.loadepi(pfile, echo, readoutFile);  % [nx*decimation ny nz nCoils]
end

imNegFreq = toppe.utils.ift3(d(:,:,:,:,1));  % -4000 Hz
imPosFreq = toppe.utils.ift3(d(:,:,:,:,2));  % +4000 Hz
bs.amp = 0.05;           % Amplitude of Fermi pulse (Gauss)
bs.freq = [-4000 4000];  % Frequency offset of Fermi pulse (Hz)
rf = toppe.readmod([datDir 'bs.mod']);
dt = 4e-6;  % RF raster time (sec)
kbs = toppe.utils.rf.calckbs(rf, bs.freq(2), dt);
pc = toppe.utils.phasecontrast3d(imPosFreq, imNegFreq);
b1 = sqrt(pc/2/kbs);   % Measured b1 in Gauss. Note factor of 2!

