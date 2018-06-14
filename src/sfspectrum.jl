import Base: +,-,*,/

"""
A container that holds a sum frequency spectrum.
"""
mutable struct SFSpectrum{T<:Number,N} <: AbstractArray{T,N}
  id::Int
  s::Array{T,N}
end

Base.size(s::SFSpectrum) = size(s.s)
Base.getindex(s::SFSpectrum, i::Int) = getindex(s.s, i)
Base.getindex(s::SFSpectrum{T,N}, I::Vararg{Int, N}) where {N,T} = getindex(s.s, I...)
Base.setindex!(s::SFSpectrum, v::Number, i::Int) = setindex!(s.s, v, i)
Base.setindex!(s::SFSpectrum{T,N}, v::Number, I::Vararg{Int, N}) where {N,T} = setindex!(s.s, v, I...)
Base.copy(s::SFSpectrum) = SFSpectrum(copy(s.id), copy(s.s))


"""
Combine all single spectra to a single one. Returns the mean counts
per second.
"""
function Base.mean!(s::SFSpectrum)
  exptime = get_attribute(s, "ccd_exposure_time")::Float64
  s.s = mean(s, 3) / exptime
  for n = ndims(s):-1:1
    if size(s, n) == 1
      s.s = squeeze(s, n)
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