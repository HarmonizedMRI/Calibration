# B1+ mapping for GE

## Bloch-Siegert B1+ mapping

Example usage:
```
sys = toppe.systemspecs('maxSlew', 10, 'gradient', 'xrm');                
N = [100 100 8];    % matrix size
FOV = [20 20 1.6];  % cm
nShots = 5;         % number of EPI segments
tr = 150;           % sequence TR (ms)
blochsiegert4ge(sys, N, FOV, nShots, tr, 'flyback', false, 'rampsamp', true);
toppe.plotseq(1,6,sys);

```

