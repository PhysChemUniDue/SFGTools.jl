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

+(s::SFSpectrum, t::SFSpectrum) = s.s + t.s
-(s::SFSpectrum, t::SFSpectrum) = s.s - t.s
*(s::SFSpectrum, a::Number)     = s.s * a
*(a::Number, s::SFSpectrum)     = a * s.s
/(s::SFSpectrum, a::Number)     = s.s / a
/(a::Number, s::SFSpectrum)     = a ./ s.s