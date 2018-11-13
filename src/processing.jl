import Statistics: mean, std



"""
`rm_events!(s::SFSpectrum, width=2)`
Remove spikes from a spectrum.
Takes a single spectrum or an array of spectra. `width` should be equal to the
width of the expected spikes (i.e. if a single spike takes up 4 pixels `width`
should be 4)

The function prints information about removed events by default. Turn off this behaviour by
setting the keyword argument `printinfo` to `false`.
"""
function rm_events!(s::SFSpectrum, width=2; kwargs...)
  s.s = rm_events!(s.s, width; kwargs...)
end

function rm_events!(s::Array{SFSpectrum}, width=2; kwargs...)
  for r in s
    r.s = rm_events!(r.s, width; kwargs...)
  end
  return s
end

function rm_events!(s::Array{Float64}, width=3; printinfo=true, minstd=5)
  r = reshape(s, (:, 1))
  dr = diff(r, dims=1)
  printinfo && (eventcounter = 0)

  threshold = minstd * std(dr)
  startidx = width + 1
  endidx = length(dr)

  for i = startidx:endidx
    lowidx = i - width

    if any(x -> x > threshold, dr[lowidx:i]) && any(x -> x < -threshold, dr[i])
      tmp = dr[lowidx:i]
      width_real = width - argmax(tmp) + 1
      r[i-width_real+1:i] .= (r[i-width_real] + r[i+1])/2
      printinfo && (eventcounter += 1)
    end
  end

  s = reshape(r, size(s))
  printinfo && println("Removed $eventcounter events.")
  return s
end


"""
`fieldcorrection!(spectrum::SFSpectrum{T,3}, width=2)`
Remove bias and dark counts from spectrum and normalize it.
Takes a single spectrum or an array of spectra. `width` should be equal to the
width of the expected spikes (i.e. if a single spike takes up 4 pixels `width`
should be 4)

The function prints information about removed events by default. Turn off this behaviour by
setting the keyword argument `printinfo` to `false`.
"""
function fieldcorrection!(spectrum::AbstractArray{T,3};
                          bias=[]::AbstractArray{T,3},
                          dark=[]::AbstractArray{T,3},
                          flat=[]::AbstractArray{T,3},
                          darkflat=[]::AbstractArray{T,3}) where T

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
    !isempty(dark) && (dark_ub = rm_offset( dark, bias))
    !isempty(flat) && (flat_ub = rm_offset( flat, bias))
    !isempty(dafl) && (dafl_ub = rm_offset( dafl, bias))
    !isempty(dark) && rm_offset!(spec,    dark_ub)
    !isempty(dafl) && rm_offset!(flat_ub, dafl_ub)

    !isempty(flat) && flatcorrection!(spec, flat_ub)

    return spec
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
function average(s::SFSpectrum{T,N}) where {T,N}
  N == 3 || error("The number of dimensions of the spectrum has to be 3.")
  exptime = get_attribute(s, "ccd_exposure_time")[1]
  if size(s, 2) == 1
    n = SFSpectrum(s.id, Array{T,1}(undef, size(s, 1)))
    n.s = mean(s, dims=3)[:,1,1] / exptime
  else
    n = SFSpectrum(s.id, Array{T,2}(undef, (size(s, 1), size(s, 2)) ))
    n.s = mean(s, dims=3)[:,:,1] / exptime
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
