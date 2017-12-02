module SFGTools

using FileIO
using DataFrames

include("sfspectrum.jl")
include("spectraio.jl")
include("processing.jl")

export  SFSpectrum,
        grab,
        list_spectra,
        load_spectra,
        get_attribute,
        get_metadata,
        rm_events!,
        rm_background!!,
        get_ir_wavelength,
        get_ir_wavenumber,
        
end
