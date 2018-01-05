using Base.Test
using DataFrames
using SFGTools

function basicfuns()
    grab(Pkg.dir("SFGTools") * "/test/sampledata/"; getall=true)
    df = list_spectra()
    id = df[:id]
    data = load_spectra(id)
    mean_counts = maximum(mean(data[5]))
    approx_wavenumber = round(mean(get_ir_wavenumber(data[1])), -1)

    (mean_counts, approx_wavenumber)
end

@test basicfuns() == (977.6, 2950.0)
