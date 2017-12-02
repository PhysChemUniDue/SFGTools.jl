"""
A container that holds a sum frequency spectrum.
"""
mutable struct SFSpectrum
  id::Int64
  s::Array{Float64}
end
