using Test
using DataFrames
using SFGTools
import Statistics: mean

const SAMPLE_DATA_DIR = joinpath(@__DIR__, "sampledata/")

function listtest()
    grab(SAMPLE_DATA_DIR; getall=true)
    grab(SAMPLE_DATA_DIR)
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
    spectrum_integer = load_spectra(63654390286607, UInt16)
    @test spectrum_integer[1][1] == 0x02be
    @test typeof(spectrum_integer) == Array{SFGTools.SFSpectrum,1}
    spectrum = load_spectra(63654390286607)
    @test typeof(spectrum) == Array{SFGTools.SFSpectrum,1}
    @test spectrum[1][1] == 702.0
    spectrum
end

function attribute_test(spectrum)
    meta = get_metadata(spectrum)
    @test typeof(meta) == Dict{String,Any}
    attr = get_attribute(spectrum, "ccd_exposure_time")
    @test attr == Any[0.5, 0.5, 0.5, 0.5, 0.5]
end

function fieldcorrection_test()
    a = Array{Float64,3}(undef, (4, 1, 2))
    a[:,1,1] = [6, 7, 8, 9]
    a[:,1,2] = [6, 8, 9, 7]
    spectrum = SFSpectrum(0, a)

    b = Array{Float64,3}(undef, (4, 1, 2))
    b[:,1,1] = [0, 1, 3, 0]
    b[:,1,2] = [2, 1, 1, 0]
    bias = SFSpectrum(0, b)

    d = Array{Float64,3}(undef, (4, 1, 2))
    d[:,1,1] = [3, 4, 3, 5]
    d[:,1,2] = [4, 2, 1, 5]
    dark = SFSpectrum(0, d)

    f = Array{Float64,3}(undef, (4, 1, 2))
    f[:,1,1] = [2, 3, 4, 3]
    f[:,1,2] = [4, 3, 6, 1]
    flat = SFSpectrum(0, f)

    l = Array{Float64,3}(undef, (4, 1, 2))
    l[:,1,1] = [1, 2, 3, 0]
    l[:,1,2] = [3, 2, 3, 2]
    darkflat = SFSpectrum(0, l)

    fieldcorrection!(spectrum, bias, dark=dark, flat=flat, darkflat=darkflat)

    @test spectrum[:,1,1] == [5.0, 8.0,  6.0, 8.0]
    @test spectrum[:,1,2] == [5.0, 10.0, 7.0, 4.0]
end

# function save_mat_test(spectra)
#     success = false
#     try
#         save_mat(tempname(), spectra)
#         success = true
#     catch
#         success = false
#     end
#     @test success == true
# end

function makespectraarray(spectrum::SFSpectrum)
    spectra = Array{SFSpectrum,1}(undef, size(spectrum,2))
    for i = 1:size(spectrum, 2)
        spectra[i] = SFSpectrum(spectrum.id, spectrum[:,i,1])
    end
    spectra
end

@testset "list_spectra Tests" begin listtest() end
@testset "load_spectra Tests" begin global spectrum = loadtest()[1] end
@testset "Attribute Tests" begin attribute_test(spectrum) end
@testset "Fieldcorrection Tests" begin fieldcorrection_test() end
# spectra = makespectraarray(spectrum)
# @testset "MAT Saving" begin save_mat_test(spectra) end
