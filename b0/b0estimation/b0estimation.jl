using MRIFieldmaps: b0init, b0map
using ROMEO: unwrap
using MAT: matread, matwrite

"""
    b0 = estimate_b0(data, echotimes; smap, outfile, b0map_kwargs...)

Estimate a (regularized and unwrapped) B0 field map from multi-echo data.

# Arguments
- `data`: One of the following:
  - `Array` of size `(dims..., nc, ne)`
    containing complex coil images (*not* k-space data)
  - `AbstractString` (e.g., `"path/to/data.mat"`) indicating a .mat file
    that contains a variable named `data`
    that is an `Array` of the form indicated above
- `echotimes`: List of echo times in seconds

# Options
- `smap = ones(dims..., nc)`: `Array` of size `(dims..., nc)`
  containing sensitivity maps for each coil;
  default is to assume uniform sensitivities
- `outfile = "b0.mat"`: `AbstractString` indicating
  where to save the estimated field map to disk;
  the resulting .mat file will contain one variable, named `b0`;
  set `outfile = ""` to skip saving the output to disk
- `b0map_kwargs...`: Additional keyword arguments passed to `MRIFieldmaps.b0map`

# Returns
- `b0`: `Array` of size `(dims...,)`
  containing the estimated unwrapped and regularized B0 field map in Hz

## Notes
- `dims...` denotes the dimensions of the data,
  i.e., `(nx, ny)` for 2D or `(nx, ny, nz)` for 3D
- `nc` is the number of coils
- `ne` is the number of echoes
"""
function estimate_b0(data, echotimes;
    smap = ones(eltype(data), size(data)[1:end-1]),
    outfile = "b0.mat",
    b0map_kwargs...
)

    # Get an initial estimate of the field map
    finit = b0init(data, echotimes; smap)

    # Compute the square root sum-of-squares coil-combined image
    # to aid with the unwrapping procedure
    imsos = sqrt.(dropdims(sum(abs2, data[:,:,:,:,1]; dims = 4); dims = 4))

    # Convert the initial field map to radians
    ΔTE = echotimes[2] - echotimes[1]
    phase = (2π * ΔTE) .* finit

    # Unwrap, then convert back to units of frequency (cycles/unit time)
    unwrapped = unwrap(phase; mag = imsos) ./ (2π * ΔTE)

    # Regularize the unwrapped field map
    (fhat,) = b0map(data, echotimes; finit = unwrapped, smap, b0map_kwargs...)

    # Save the estimated field map to file
    outfile == "" || matwrite(outfile, Dict("b0" => fhat))

    return fhat

end

function estimate_b0(datafile::AbstractString, args...; kwargs...)

    data = matread(datafile)["data"]
    return estimate_b0(data, args...; kwargs...)

end
