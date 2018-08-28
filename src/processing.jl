import Statistics: mean, std

const VIS_WAVELENGTH = 512.4  # According to service protocol of May 2018


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
Combine all single spectra to a single one. Returns the mean counts
per second.
"""
function average(s::SFSpectrum{T,N}) where {T,N}
  N == 3 || error("The number of dimensions of the spectrum has to be 3.")
  exptime = get_attribute(s, "ccd_exposure_time")::Float64
  if size(s, 2) == 1
    n = SFSpectrum(s.id, Array{T,1}(size(s, 1)))
    n.s = mean(s.s, 3)[:,1,1] / exptime
  else
    n = SFSpectrum(s.id, Array{T,2}(size(s, 1, 2)))
    n.s = mean(s.s, 3)[:,:,1] / exptime
  end
  n
end


"""
Substract the blind counts stored in `bg` from each of the spectra of `s`.
"""
function rm_blindcounts!(s::SFSpectrum, bg::SFSpectrum)
    # Check if size is correct
    spectrum_size = (size(s, 1), size(s, 2))
    bg_size = (size(bg, 1), size(bg, 2))
    spectrum_size == bg_size || error("Sizes of spectrum and background do not match. The spectrum is of size $spectrum_size whereas the background is of size $(bg_size).")

    # Mean counts of background
    bg_counts = mean(bg, dims=3)

    # Subtract the background from each spectrum of s
    for i = 1:size(s, 3)
        s.s[:,:,i] -= bg_counts
    end
    s
end


"""
Remove Background from the first spectrum passed to the function. →SFSpectrum
Takes a spectrum or an array of spectra `s` and a background spectrum `bg`.
"""
function rm_background!(s::SFSpectrum, bg::SFSpectrum)
    # Make a copy of the background spectrum because we don't want to
    # change it inside this function.
    bgc = copy(bg)
    a = [s, bgc]
    for i = 1:2
        if ndims(a[i]) > 3
            @show size(a[i])
            println("Expected ndims(s) and ndims(bg) to be 1, 2 or 3.")
            return
        elseif ndims(a[i]) == 3
            mean!(a[i])
        end
    end
    s.s -= bgc.s
end

function rm_background!(a::Array{SFSpectrum}, bg::SFSpectrum)
    for i = 1:length(a)
        a[i].s = rm_background!(a[i], bg)
    end
    return a
end


"""
Return the wavelength of the detected light in nm.

The appropriate calibration curve of the spectrometer is automatically selected via the timestamp of the
spectrum. If you want to select a specific calibration set the `date` keyword argument to something like
`date=DateTime("2011-11-11")` to select the calibration curve that was valid on that specific date.
"""
function get_wavelength(s::SFSpectrum;
                        date = DateTime(get_attribute(s, "timestamp")))

    λ0 = get_attribute(s, "spectrometer_wavelength")
    x_binning = get_attribute(s, "x_binning")
    num_points = N_PIXEL / x_binning

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

    wl_range = linspace(first_pixel_wl, last_pixel_wl, num_points)
    wl = collect(wl_range)

end

"""
Return the wavenumbers in inverse centimeters

The appropriate calibration curve of the spectrometer is automatically selected via the timestamp of the
spectrum. If you want to select a specific calibration set the `date` keyword argument to something like
`date=DateTime("2011-11-11")` to select the calibration curve that was valid on that specific date.

You can change the default wavelength for the visible light by passing a value to the `vis` keyword argument.
"""
function get_wavenumber(s::SFSpectrum; date=DateTime(get_attribute(s, "timestamp")))
    1 ./ get_wavelength(s; kwargs...) * 1e7
end

"""
Return the wavelength of the corresponding infrared light in nm.
"""
function get_ir_wavelength(s::SFSpectrum; vis=VIS_WAVELENGTH, date=DateTime(get_attribute(s, "timestamp")))

    function sf2ir(sf, vis)
        1 ./ (1 ./ sf - 1 ./ vis)
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
Get the parameter that changes during the measurement as a string
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
    p = get_attribute(s, "pump_dl_position")
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

"""
Automated Processing.
"""
function quick_process(d::Array{SFSpectrum}; width=5)
    # 1. Get the variable of the Dataset
    vars = get_variables(d)

    # 2. Processing
    rm_events!.(d, width)
    mean!.(d)

    λ = get_ir_wavelength(d[1])

    z = zeros(size(d[1].s, 1), length(d), size(d[1].s, 2))
    for i = 1:length(d)
        z[:,i,:] = d[i].s
    end

    return λ, vars, z
end
