import Base: mean, mean!, size, ndims
import Base: +,-,*,/
import Base: copy, show

"""
A container that holds a sum frequency spectrum.
"""
mutable struct SFSpectrum
  id::Int64
  s::Array{Float64}
end

size(s::SFSpectrum) = size(s.s)
ndims(s::SFSpectrum) = length(size(s))

"""
Combine all single spectra to a single one. Returns the mean counts
per second.
"""
function mean(s::SFSpectrum)
    mean_spec = copy(s)
    mean!(mean_spec)
    return mean_spec
end

function mean!(s::SFSpectrum)
  exptime = get_attribute(s, "ccd_exposure_time")::Float64
  s.s = mean(s.s, 3) / exptime
  for n = ndims(s.s):-1:1
    if size(s.s, n) == 1
      s.s = squeeze(s.s, n)
    end
  end
  return s
end

+(s::SFSpectrum, t::SFSpectrum) = s.s + t.s
-(s::SFSpectrum, t::SFSpectrum) = s.s - t.s
*(s::SFSpectrum, a::Number)     = s.s * a
*(a::Number, s::SFSpectrum)     = a * s.s
/(s::SFSpectrum, a::Number)     = s.s / a
/(a::Number, s::SFSpectrum)     = a ./ s.s


copy(s::SFSpectrum) = SFSpectrum(copy(s.id), copy(s.s))

show(io::IO, d::SFSpectrum) = print(io,"$(typeof(d))($(d.id), $(typeof(d.s)) size $(size(d.s)))")