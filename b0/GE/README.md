# B0 mapping for GE

## 3D FLASH

Example usage:
```
sys = toppe.systemspecs('maxSlew', 15, 'gradient', 'xrm');                
N = [100 100 8];    % matrix size
FOV = [20 20 1.6];  % cm
flip = 5;           % degrees
DTE = [0 2.3];      % (ms) acquire 2 images with TE extended by 0 and 2.3 ms
b04ge(sys, N, FOV, flip, DTE);
toppe.plotseq(1,4,sys);

```

