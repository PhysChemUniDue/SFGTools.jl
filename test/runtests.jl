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
    @test typeof(spectrum) == Array{SFGTools.SFSpectrum,1}
    spectrum
end

function attribute_test(spectrum)
    meta = get_metadata(spectrum)
    @test typeof(meta) == Dict{String,Any}
    attr = get_attribute(spectrum, "ccd_exposure_time")
    @test attr == 0.5
end

function blindcounts_test()
    blind_spectrum = SFSpectrum(123, rand(512,1,1000) .+ 700.0)
    spectrum = SFSpectrum(456, rand(512,1,5) .+ 710.0)
    rm_blindcounts!(spectrum, blind_spectrum)
    @test 9.0 < mean(spectrum.s) < 11.0
end

function save_mat_test(spectra)
    success = false
    try
        save_mat("/tmp/" * randstring(), spectra)
        success = true
    catch
        success = false
    end
    @test success == true
end

@testset "list_spectra Tests" begin listtest() end
@testset "load_spectra Tests" begin global spectrum = loadtest() end
@testset "Attribute Tests" begin attribute_test(spectrum) end
@testset "Blindcounts Removal" begin blindcounts_test() end
spectra = makespectraarray(spectrum[1])
@testset "MAT Saving" begin save_mat_test(spectra) end

function makespectraarray(spectrum::SFSpectrum)
    spectra = Array{SFSpectrum,1}(size(spectrum,2))
    for i = 1:size(spectrum, 2)
        spectra[i] = SFSpectrum(spectrum.id, spectrum[:,i,1])
    end
    spectra
end

