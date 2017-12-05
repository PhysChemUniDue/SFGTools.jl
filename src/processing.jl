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


function get_ir_wavelength(s::SFSpectrum)
    const N_PIXEL = 512
    const VIS_WAVELENGTH = 512.8
    const RANGE = 210

    function sf2ir(sf, vis)
        1 / (1/sf - 1/vis)
    end

    sf_wavelength = get_attribute(s, "spectrometer_wavelength")
    x_binning = get_attribute(s, "x_binning")
    num_points = N_PIXEL/x_binning
    ir_central = sf2ir(sf_wavelength, VIS_WAVELENGTH)
    ir_wavelength = linspace(ir_central - RANGE, ir_central + RANGE, num_points)
end


function get_ir_wavenumber(s::SFSpectrum)
    1 / get_ir_wavelength(s) * 1e7
end
