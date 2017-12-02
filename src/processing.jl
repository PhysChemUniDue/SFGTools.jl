"""
Remove cosmic events from SFSpectrum
"""
function rm_events!(s::SFSpectrum)
  s.s = rm_events!(s.s)
end

# function rm_events!(s::Array{SFSpectrum})
#   for r in s
#     r.s = rm_events!(r.s)
#   end
#   return s
# end

function rm_events!(s::AbstractArray)
  r = reshape(s, (:, 1))
  dr = diff(r)
  threshold = 4 * std(dr)
  for i = 3:length(dr)-1
    if dr[i-1] > threshold && dr[i] < -threshold
      println(i)
      r[i] = (r[i-2] + r[i+2])/2
    end
  end
  s = reshape(r, size(s))
  # figure()
  # plot(dr)
  # hlines([threshold, -threshold], 1, length(dr))
  return s
end


"""
Remove Background →SFSpectrum
"""
function rm_background!(s::SFSpectrum, bg::SFSpectrum)
    # Make a copy of the background spectrum because we don't want to
    # change it inside this function.
    bgc = copy(bg)
    a = [s, bgc]
    for i = 1:2
        if ndims(a[i]) != 2 && ndims(a[i]) != 3
            @show i
            @show size(a[i])
            println("Expected ndims(s) and ndims(bg) to be 2 or 3.")
            return
        elseif ndims(a[i]) == 3
            mean!(a[i])
        end
    end
    s.s -= bgc.s
end

function rm_background!(a::Array{SFSpectrum}, bg::SFSpectrum)
    for i = 1:length(a)
        rm_background!(a[i], bg)
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
