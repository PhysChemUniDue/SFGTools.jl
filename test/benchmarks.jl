using BenchmarkTools
using SFGTools

println("Listing Spectra:")
df = @btime list_spectra(inexact="twinscan", group=true)

println("Loading 7 Spectra (tiff):")
data = @btime load_spectra($df[:id][1])
println("Sizes: $(size(data[1]))\n")

println("Loading single Spectrum (tiff)")
data = @btime load_spectra($df[:id][1][1])
println("Size: $(size(data))\n")

println("Event Removal:")
@btime rm_events!.($data)

println("Get Metadata:")
@btime get_metadata($data[1])
