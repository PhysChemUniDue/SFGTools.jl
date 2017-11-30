"""
Remove cosmic events from SFSpectrum
"""
function remove_events!(s::SFSpectrum)
  s.spectrum = remove_events!(s.spectrum)
end

function remove_events!(s::Array{SFSpectrum})
  for r in s
    r.spectrum = remove_events!(r.spectrum)
  end
  return s
end

function remove_events!(s::AbstractArray)
  # Make the derivative of the
  r = reshape(s, (:, 1))
  dr = diff(r)
  threshold = 3 * std(dr)
  for i = 2:length(dr)
    if dr[i-1] > threshold && dr[i] < -threshold
      println(i)
      r[i] = (r[i-1] + r[i+1])/2
    end
  end
  s = reshape(r, size(s))
  return s
end
