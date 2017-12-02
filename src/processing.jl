"""
Remove cosmic events from SFSpectrum
"""
function remove_events!(s::SFSpectrum)
  s.s = remove_events!(s.s)
end

# function remove_events!(s::Array{SFSpectrum})
#   for r in s
#     r.s = remove_events!(r.s)
#   end
#   return s
# end

function remove_events!(s::AbstractArray)
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
function rmbackground!(s::SFSpectrum, bg::SFSpectrum)
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

function rmbackground!(a::Array{SFSpectrum}, bg::SFSpectrum)
    for i = 1:length(a)
        rmbackground!(a[i], bg)
    end
    return a
end
