__precompile__()
module SFGTools

using DataFrames
# using MAT

include("sfspectrum.jl")
include("spectraio.jl")
include("processing.jl")

const N_PIXEL = 512

export  SFSpectrum,
        average,
        grab,
        list_spectra,
        load_spectra,
        get_attribute,
        get_metadata,
        is_attribute,
        rm_events!,
        rm_background!,
        rm_blindcounts!,
        get_wavelength,
        get_wavenumber,
        get_ir_wavelength,
        get_ir_wavenumber,
        get_pump_delay,
        dlpos2t,
        get_variables,
        quick_process
        # save_mat

end
