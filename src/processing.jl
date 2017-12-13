using PyPlot

"""
Remove spikes from a spectrum.
Takes a single spectrum or an array of spectra. `width` should be equal to the
width of the expected spikes (i.e. if a single spike takes up 4 pixels `width`
should be 4)
"""
function rm_events!(s::SFSpectrum, width=2)
  s.s = rm_events!(s.s, width)
end

function rm_events!(s::Array{SFSpectrum}, width=2)
  for r in s
    r.s = rm_events!(r.s, width)
  end
  return s
end

function rm_events!(s::Array{Float64}, width=2)
  r = reshape(s, (:, 1))
  dr = diff(r)

  threshold = 4 * std(dr)
  startidx = width + 2
  endidx = length(dr) - width

  for i = startidx:endidx
    lowidx = i - width
    uppidx = i + width - 1

    if any(x -> x > threshold, dr[lowidx:i]) && any(x -> x < -threshold, dr[i:uppidx])
      r[i] = (r[lowidx-1] + r[uppidx+1])/2
    end
  end

  s = reshape(r, size(s))
  return s
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


function get_wavelength(s::SFSpectrum)
    const N_PIXEL = 512

    λ0 = get_attribute(s, "spectrometer_wavelength")
    x_binning = get_attribute(s, "x_binning")
    num_points = N_PIXEL / x_binning

    # Calibration Parameters were determined in 2017-12-13_SpectrometerCalibration.ipynb
    # Pixel offset
    Δp = -0.000338165 * λ0^2 + 0.250245 * λ0 - 32.0725
    # Wavelenths per pixel
    pλ = -3.36321e-8 * λ0^2 + 4.9358e-6 * λ0 + 0.0192305

    # Central Wavelength
    λc = λ0 - Δp * pλ

    # First Pixel Wavelength:
    first_pixel_wl = λc - (num_points/2 - 0.5) * pλ
    last_pixel_wl  = λc + (num_points/2 - 0.5) * pλ

    wl_range = linspace(first_pixel_wl, last_pixel_wl, num_points)
    wl = collect(wl_range)

end


function get_ir_wavelength(s::SFSpectrum; vis=512.6)

    function sf2ir(sf, vis)
        1 ./ (1./sf - 1/vis)
    end

    sf_wavelength = get_wavelength(s)
    ir_wavelength = sf2ir(sf_wavelength, vis)
end


function get_ir_wavenumber(s::SFSpectrum)
    1 ./ get_ir_wavelength(s) * 1e7
end
