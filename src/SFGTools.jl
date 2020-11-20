__precompile__()
module SFGTools

using DataFrames
# using MAT

include("sfspectrum.jl")
include("spectraio.jl")
include("processing.jl")

const N_PIXEL = 512
const VIS_WAVELENGTH = 512.4  # According to service protocol of May 2018

export  SFSpectrum,
        average,
        grab,
        list_spectra,
        load_spectra,
        fieldcorrection!,
        get_attribute,
        get_metadata,
        is_attribute,
        rm_events!,
        get_wavelength,
        get_wavenumber,
        get_ir_wavelength,
        get_ir_wavenumber,
        get_pump_delay,
        dlpos2t,
        get_variables,
        pixelshift,
        pixelshift!,
        findpixelshift,
        savejson,
        # save_mat,
        save_data

end
