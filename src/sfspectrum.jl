"""
A container that holds a sum frequency spectrum.
"""
mutable struct SFSpectrum
  id::Int64
  name::AbstractString
  dir::AbstractString
  spectrum::Array{Float64}
end
