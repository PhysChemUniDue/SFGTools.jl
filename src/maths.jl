import Base: mean, mean!, size, ndims
import Base: +,-,*,/
import Base: copy

size(s::SFSpectrum) = size(s.s)
ndims(s::SFSpectrum) = length(size(s))

"""
Combine all single spectra to a single one. Returns the mean counts
per second.
"""
function mean(s::SFSpectrum)
    exptime = get_attribute(s, "ccd_exposure_time")
    squeeze(mean(s.s, 3), 3) / exptime
end

function mean!(s::SFSpectrum)
  s.s = mean(s)
end

+(s::SFSpectrum, t::SFSpectrum) = s.s + t.s
-(s::SFSpectrum, t::SFSpectrum) = s.s - t.s
*(s::SFSpectrum, a::Number)     = s.s * a
*(a::Number, s::SFSpectrum)     = a * s.s
/(s::SFSpectrum, a::Number)     = s.s / a
/(a::Number, s::SFSpectrum)     = a ./ s.s


copy(s::SFSpectrum) = SFSpectrum(copy(s.id), copy(s.s))
