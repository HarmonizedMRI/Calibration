# Unwrapped, Regularized B0 Field Map Estimation

[b0estimation.jl](b0estimation.jl) defines the function `estimate_b0`
for estimating an unwrapped, regularized B0 field map.


## Getting Started

`estimate_b0` depends on three Julia packages:
[MRIFieldmaps.jl](https://github.com/MagneticResonanceImaging/MRIFieldmaps.jl),
[ROMEO.jl](https://github.com/korbinian90/ROMEO.jl), and
[MAT.jl](https://github.com/JuliaIO/MAT.jl).
These packages can be installed
by doing the following from the Julia REPL:
```julia-repl
julia> ] # enter Julia's package prompt
pkg> add MRIFieldmaps ROMEO MAT
```


## Usage

For each new Julia session,
you must first `include` [b0estimation.jl](b0estimation.jl)
to load the necessary packages
and to define the function `estimate_b0`:
```julia
include("path/to/Calibration/b0/b0estimation/b0estimation.jl")
```
This needs to be done just once per Julia session.
Then the following will estimate a field map
given appropriate inputs:
```julia
b0 = estimate_b0(data, echotimes; smap, outfile)
```
See the documentation for a description of the inputs.


## Documentation

See the docstring in [b0estimation.jl](b0estimation.jl),
or run the following from the Julia REPL
after `include`ing [b0estimation.jl](b0estimation.jl):
```julia-repl
julia> ? # enter Julia's help mode
help?> estimate_b0
```
