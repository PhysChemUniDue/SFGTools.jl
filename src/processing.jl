import Statistics: mean, std
import Dierckx: Spline1D, evaluate
import Optim: BFGS, optimize


"""
`rm_events!(s::SFSpectrum; width=2, minstd=5)`
Remove spikes from a spectrum and return the number of events removed.

Takes a single spectrum or an array of spectra. `width` should be equal to the
width of the expected spikes (i.e. if a single spike takes up 4 pixels `width`
should be 4). `minstd` denotes the factor of which the peak has to be above
the standard deviation of the spectrum counts to be recognized as an event.
"""
function rm_events!(s::SFSpectrum; width=3, minstd=5)
  num_removed_events = rm_events!(s.s; width=width, minstd=minstd)
end

function rm_events!(s::Array{Float64}; width=3, minstd=5)

# Behaves like `diff` but is twice as fast
function mydiff(a)
    b = typeof(a)(undef, length(a)-1)
        for i = 1:length(a)-1
           b[i] = a[i+1] - a[i]
        end
    b
end

  r = reshape(s, :)
  dr = mydiff(r)
  num_removed_events = 0

  threshold = minstd * std(dr)
  startidx = width + 1
  endidx = length(dr)
  lowidx = 0

  for i = startidx:endidx
    lowidx = i - width

    if any(dr[j] > threshold for j in lowidx:i) && dr[i] < -threshold
      tmp = dr[lowidx:i]
      width_real = width - argmax(tmp) + 1
      r[i-width_real+1:i] .= (r[i-width_real] + r[i+1])/2
      num_removed_events += 1
    end
  end

  s .= reshape(r, size(s))

  num_removed_events
end


"""
`pixelshift!(a::AbstractArray{T,N}, dim::Int, shift::Real, amplification::Real=1.0)`
Shift the array a the number of pixels `shift` alongside dimension `dim` and
amplify the values by `amplification`.

Example:
```
julia> a = [2.0, 4.0, 6.0, 2.0];
julia> pixelshift!(a, 1, 0.5, 1.0)
4-element Array{Float64,1}:
 3.0
 5.0
 4.0
 2.0
```
"""
function pixelshift!(a::AbstractArray{T,N}, dim::Int, shift::Real, amplification::Real=1.0) where {T <: Real, N}

    N > 2 && error("pixelshift not defined for the case of $N dimensions")

    pixels = 1:size(a, dim)
    pixels_shifted = pixels .+ shift

    if N == 1
        spl = Spline1D(pixels, a, k=1)
        a .= evaluate(spl, pixels_shifted) .* amplification
    elseif N >= 2
        for i = 1:size(a', dim)
            dim == 1 && (spl = Spline1D(pixels, a[:,i], k=1))
            dim == 1 && (a[:,i] .= evaluate(spl, pixels_shifted) .* amplification)
            dim == 2 && (spl = Spline1D(pixels, a[i,:], k=1))
            dim == 2 && (a[i,:] .= evaluate(spl, pixels_shifted) .* amplification)
        end
    end

    a

end

function pixelshift(a::AbstractArray{T,N}, dim::Int, shift::Real, amplification::Real=1.0) where {T <: Real, N}

    b = deepcopy(a)
    pixelshift!(b, dim, shift, amplification)

end


"""
`findpixelshift(a::AbstractArray{T}, b::AbstractArray{T}, dim::Int)`
Find the optimal `shift` and `amplification` values for the `pixelshift`
function.

Array `b` is shifted alongside `a` in dimension `dim` to find the optimal
shift parameters. The result is on `Optim.jl` result. The shift can be accessed
by `result.minimizer[1]` and the amplitude by `result.minimizer[2]`.
"""
function findpixelshift(a::AbstractArray{T}, b::AbstractArray{T}, dim::Int) where T <: Number
    minfun(x) = mean((a .- pixelshift(b, dim, x[1], x[2])).^2)
    result = optimize(minfun, [0.0, 1.0], BFGS())
end


"""
`fieldcorrection!(spectrum::AbstractArray{T,3}; bias=[]::AbstractArray{T,3}, dark=[]::AbstractArray{T,3}, flat=[]::AbstractArray{T,3}, darkflat=[]::AbstractArray{T,3})`
Remove bias and dark counts from spectrum and normalize it with a flat spectrum.

All spectra provided have to be 3D (including bias etc.).

The algorithm works as follows where ``s`` is the spectrum ``b`` the bias,
``d`` the dark spectrum, ``f`` the flat spectrum and ``l`` the dark flat
spectrum. The mean indicated by a ``< >`` is taken along the third dimension
of the spectra.

``s_b = s - <b>``

``d_b = d - <b>``

``f_b = f - <b>``

``l_b = l - <b>``

``s_{b,d} = s_b - <d_b>``

``f_{b,d} - f_b - <l_b>``

``f_{b,d,norm} = <f_{b,d}> / <f_{b,d}>_{max}``

``s_{b,d,f} = s_{b,d} / f_{b,d,norm}``

If no bias is provided ``<b>`` equals zero.
"""
function fieldcorrection!(spectrum::AbstractArray{T,3};
                          bias=Array{Float64}(undef,0,0,0)::AbstractArray{T,3},
                          dark=Array{Float64}(undef,0,0,0)::AbstractArray{T,3},
                          flat=Array{Float64}(undef,0,0,0)::AbstractArray{T,3},
                          darkflat=Array{Float64}(undef,0,0,0)::AbstractArray{T,3}) where T

    function rm_offset!(s, offset)
        offset_m = mean(offset, dims=3)
        for i = 1:size(s, 3)
            s[:,:,i] .-= offset_m[:,:,1]
        end
        s
    end

    function rm_offset(s, offset)
        offset_m = mean(offset, dims=3)
        n = deepcopy(s)
        for i = 1:size(s, 3)
            n[:,:,i] .-= offset_m[:,:,1]
        end
        n
    end

    function flatcorrection!(s, flat)
        flatmean = mean(flat, dims=3)
        flatmean ./= maximum(flatmean)
        for i = 1:size(s, 3)
            s[:,:,i] ./= flatmean[:,:,1]
        end
        s
    end

    # For readability make same length as all
    dafl = darkflat
    spec = spectrum

    !isempty(bias) &&            rm_offset!(spec, bias)

    if !isempty(bias)
        !isempty(dark) && (dark_ub = rm_offset( dark, bias))
        !isempty(flat) && (flat_ub = rm_offset( flat, bias))
        !isempty(dafl) && (dafl_ub = rm_offset( dafl, bias))
    elseif isempty(bias) && !isempty(dark)
        !isempty(dark) && (dark_ub = dark)
        !isempty(flat) && (flat_ub = flat)
        !isempty(dafl) && (dafl_ub = dafl)
    end

    !isempty(dark) && rm_offset!(spec,    dark_ub)
    !isempty(dafl) && rm_offset!(flat_ub, dafl_ub)

    !isempty(flat) && flatcorrection!(spec, flat_ub)

    return spec
end


"""
    get_exposure_time(s)
    
    Returns exposure time for the spectrum s.

"""
function get_exposure_time(s::SFSpectrum{T,N}) where {T,N}
    get_attribute(s, "ccd_exposure_time")[1]
end


"""
    average(s)

Combine all spectra of a series in `s` to a single one. Returns the mean counts
per second. In more detail we take the mean of the third dimension and
divide that by the exposure time.

Note that there is also a bias of around
710 count for each readout that has to be subtracted first to get the true
mean counts per second.

# Example
```julia-repl
julia> s
3×2×2 SFSpectrum{Float64,3}:
[:, :, 1] =
 0.6  1.0
 0.2  0.2
 0.2  0.2

[:, :, 2] =
 0.6  0.6
 1.0  0.8
 0.2  0.4

julia> average(s)
3×2 SFSpectrum{Float64,2}:
 0.6  0.8
 0.6  0.5
 0.2  0.3

julia> get_attribute(s, "ccd_exposure_time")[1]
1.0
```
"""
function average(s::SFSpectrum{T,N}; combine = true) where {T,N}
  N == 3 || error("The number of dimensions of the spectrum has to be 3.")
  exptime = get_attribute(s, "ccd_exposure_time")[1]
  if size(s, 2) == 1
      if combine == true
          n = SFSpectrum(s.id, Array{T,1}(undef, size(s, 1)))
          n.s = mean(s, dims=3)[:,1,1] / exptime
      else
          n = SFSpectrum(s.id, s.s)
          n.s = s.s ./ exptime
      end
  else
      if combine == true
          n = SFSpectrum(s.id, Array{T,2}(undef, (size(s, 1), size(s, 2)) ))
          n.s = mean(s, dims=3)[:,:,1] / exptime
      else
          n = SFSpectrum(s.id, s.s)
          n.s = s.s ./ exptime
      end
  end
  n
end


"""
Return the wavelength of the detected light in nm.

The appropriate calibration curve of the spectrometer is automatically selected via the timestamp of the
spectrum. If you want to select a specific calibration set the `date` keyword argument to something like
`date=DateTime("2011-11-11")` to select the calibration curve that was valid on that specific date.
"""
function get_wavelength(s::SFSpectrum;
                        date = DateTime(get_attribute(s, "timestamp")[1]))

    λ0 = get_attribute(s, "spectrometer_wavelength")[1]
    x_binning = get_attribute(s, "x_binning")[1]
    num_points = N_PIXEL ÷ x_binning

    if date < DateTime("2018-05-14")
        # Calibration Parameters were determined in 2017-12-13_SpectrometerCalibration.ipynb
        # Pixel offset
        Δp = -0.000338165 * λ0^2 + 0.250245 * λ0 - 32.0725
        # Wavelenths per pixel
        pλ = -3.36321e-8 * λ0^2 + 4.9358e-6 * λ0 + 0.0192305
    else
        # Calibration parameters were determined on 2018-05-14T17:35:18.995
        # Pixel offset
        Δp = -0.0002511567068034343 * λ0^2 + 0.17459249760169204 * λ0 - -12.204165364873228
        # Wavelenths per pixel
        pλ = -2.8767121934619073e-8 * λ0^2 + 3.532087562855135e-7 * λ0 + 0.02024248616034501
    end

    # Central Wavelength
    λc = λ0 - Δp * pλ

    # First Pixel Wavelength:
    first_pixel_wl = λc - (num_points/2 - 1/2) * pλ * x_binning
    last_pixel_wl  = λc + (num_points/2 - 1/2) * pλ * x_binning

    wl_range = range(first_pixel_wl, stop=last_pixel_wl, length=num_points)
    wl = collect(wl_range)

end

"""
Return the wavenumbers in inverse centimeters

The appropriate calibration curve of the spectrometer is automatically selected via the timestamp of the
spectrum. If you want to select a specific calibration set the `date` keyword argument to something like
`date=DateTime("2011-11-11")` to select the calibration curve that was valid on that specific date.

You can change the default wavelength for the visible light by passing a value to the `vis` keyword argument.
"""
function get_wavenumber(s::SFSpectrum; date=DateTime(get_attribute(s, "timestamp")[1]))
    1 ./ get_wavelength(s; date=date) * 1e7
end

"""
Return the wavelength of the corresponding infrared light in nm.
"""
function get_ir_wavelength(s::SFSpectrum; vis=VIS_WAVELENGTH, date=DateTime(get_attribute(s, "timestamp")[1]))

    function sf2ir(sf, vis)
        1 ./ (1 ./ sf .- 1 ./ vis)
    end

    sf_wavelength = get_wavelength(s; date=date)
    ir_wavelength = sf2ir(sf_wavelength, vis)
end

"""
Return the wavenumber of the corresponding infrared light in inverse centimeters.

The appropriate calibration curve of the spectrometer is automatically selected via the timestamp of the
spectrum. If you want to select a specific calibration set the `date` keyword argument to something like
`date=DateTime("2011-11-11")` to select the calibration curve that was valid on that specific date.

You can change the default wavelength for the visible light by passing a value to the `vis` keyword argument.
"""
function get_ir_wavenumber(s::SFSpectrum; kwargs...)#vis=VIS_WAVELENGTH, date=DateTime(get_attribute(s, "timestamp")))
    1 ./ get_ir_wavelength(s; kwargs...) * 1e7
end


"""
Get the parameters where there is change during the measurement as a dict
"""
function get_variables(d::Array{SFSpectrum})

    dict = Dict{String, AbstractArray}()

    # Check if there is data
    length(d) == 0 && return dict

    keylist = d[1] |> get_metadata |> keys

    for k in keylist
        vals = get_attribute(d, k)
        uvals = unique(vals)
        if length(uvals) > 1
            dict[k] = uvals
        end
    end

    return dict

end


"""
Get the pump delay time in ps for a given SFSpectrum
"""
function get_pump_delay(s::SFSpectrum)
    p = get_attribute(s, "pump_dl_position")[1]
    dlpos2t(p)
end

"""
Return the delay time for the infrared pump in ps
"""
function dlpos2t(p)
    c = 299792458.0
    l = 0.61
    n = 125000.0
    p * l / (n * c) * 1e12
end
