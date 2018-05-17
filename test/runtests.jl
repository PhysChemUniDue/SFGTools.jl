using Base.Test
using DataFrames
using SFGTools

const SAMPLE_DATA_DIR = Pkg.dir("SFGTools") * "/test/sampledata/"

function listtest()
    grab(SAMPLE_DATA_DIR; getall=true)
    df = list_spectra()
    @test size(df,1) == 118
    df = list_spectra(date=(2018,2,15))
    @test size(df,1) == 118
    df = list_spectra(inexact="HDT")
    @test size(df,1) == 78
    df = list_spectra(exact="HDT_180109A_EksplaScan_FR")
    @test size(df,1) == 25
    df = list_spectra(group=true)
    @test size(df,1) == 7
end

function loadtest()
    spectrum = load_spectra(63654390286607)
    @test typeof(spectrum) == SFSpectrum
    spectrum
end

function attribute_test(spectrum)
    meta = get_metadata(spectrum)
    @test typeof(meta) == Dict{String,Any}
    attr = get_attribute(spectrum, "ccd_exposure_time")
    @test attr = 0.5
end

@testset "list_spectra Tests" begin listtest() end
@testset "load_spectra Tests" begin spectrum = loadtest() end
@testset "Attribute Tests" begin attribute_test(spectrum) end
