import Base: mean, mean!, size
import Base: +,-,*,/

mean(s::SFSpectrum, dim=3) = squeeze(mean(s.spectrum, dim), dim)
size(s::SFSpectrum)        = size(s.spectrum)

function mean!(s::SFSpectrum)
  s.spectrum = mean(s, dim)
end

+(s::SFSpectrum, t::SFSpectrum) = s.spectrum + t.spectrum
-(s::SFSpectrum, t::SFSpectrum) = s.spectrum - t.spectrum
*(s::SFSpectrum, a::Number)     = s.spectrum * a
*(a::Number, s::SFSpectrum)     = a * s.spectrum
/(s::SFSpectrum, a::Number)     = s.spectrum / a
/(a::Number, s::SFSpectrum)     = a ./ s.spectrum
