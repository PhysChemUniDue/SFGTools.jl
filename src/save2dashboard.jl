using Dates,HDF5 
""" 
Convert a date (yyyy,mm,dd) of type Tuple{Int64,Int64,Int64} to a yyyymmdd string.
"""
function date2dashboard(date::Tuple{Int64,Int64,Int64})
    if ndigits(date[1]) == 4
        yyyy = string(date[1])
    else 
        error("$string(date[1]) is not a year. Format must be yyyy.")
    end
    
    if ndigits(date[2]) == 2
        mm = string(date[2])
    elseif  ndigits(date[2]) == 1
        mm = "0"*string(date[2])
    else 
        error("$string(date[2]) is not a month. Format must be either mm or m.")
    end
    
    if ndigits(date[3]) == 2
        dd = string(date[3])
    elseif  ndigits(date[3]) == 1
        dd = "0"*string(date[3])
    else 
        error("$string(date[3]) is not a day. Format must be either dd or d.")
    end
    
    return yyyy*mm*dd
end

"""
Change the - or + to ⁻ or ⁺
"""

function convert2unicode(scantype)
    if occursin("-",scantype) 
        return replace(scantype,"-"=> "⁻")
    elseif occursin("+",scantype)
        return replace(scantype,"+" => "⁺")
    else 
        return scantype
    end
end
        

""" Save the data in a HDF5 File. \n

Register a new Sample name or choose a sample name from the SFGDashboard Sample Dropdown, e.g. "CaAra" or "CaAra d4". If you are not sure visit --> 132.252.80.59:8050/\n
Example:\n
    \tsample = "CaAra" \n
    \tsample_prep = "20200504"\n
    \tfoldername =  "CaAra_20200504_1_DL-Scan_r--pumped_ppp"\n
    save_data(sample,32,"ppp","DL d-")\n
    save_data(  sample::String, surface_density_value::Int, polarisation_comb::String, scan::String;\n
                date=date, sample_prep= sample_prep, foldername= foldername, raw_spectra = raw,\n
                sigmatrix = sig03, refmatrix = ref03, probe_wavenumbers = ν, delay_time = dltime_sorted,\n
                pump_wavenumbers = nosthing, mode_name = mode_name, \n
                sig_bleaches = [mean(sig03[:,pixel[i]], dims=2)[:,1] for i in 1:length(mode_name)],\n
                ref_bleaches = [mean(ref03[:,pixel[i]], dims=2)[:,1] for i in 1:length(mode_name)],\n
                add_comment= "", \n 
                save_path = "./"
    )
"""
function save_data( sample::String, surface_density_value::Int, polarisation_comb::String, scan::String; 
                    date=Main.date, sample_prep= Main.sample_prep, foldername= Main.foldername, raw_spectra = Main.raw, 
                    sigmatrix = Main.sig03, refmatrix = Main.ref03, probe_wavenumbers = Main.ν, delay_time = Main.dltime_sorted,
                    pump_wavenumbers = nothing, mode_name = Main.mode_name, 
                    sig_bleaches = [mean(Main.sig03[:,Main.pixel[i]], dims=2)[:,1] for i in 1:length(Main.mode_name)],
                    ref_bleaches = [mean(Main.ref03[:,Main.pixel[i]], dims=2)[:,1] for i in 1:length(Main.mode_name)],
                    add_comment= "",
                    save_path = "./" 
        )

    # Check if date has the right type
    if typeof(date) == String
        directory = date * "/" * foldername
        dashboard_date = date
    elseif  typeof(date) !== String 
        dashboard_date = date2dashboard(date)
        directory = dashboard_date * "/" * foldername
    elseif length(date) !== 8
        error("$date has not the type  yyyymmdd")
    end

    # sample surface density
    surface_density = "$surface_density_value mNm⁻¹"
        
    # scan type (delay or wavenumber scan)
    if occursin("dl",lowercase(scan)) == true
        scan_type = "delay_scan"
        pump_resonance= convert2unicode(split(scan)[2])*"pumped"
    elseif occursin("wl",lowercase(scan)) == true || occursin("wn",lowercase(scan)) == true
        scan_type = "wavenumber_scan"
    else 
        error("""Scan type cannot be $scan !!\n
        Example for delay scan with d- pumped scan ="dl d-"\n
        Example for wavelength scan scan = "wl".
        """)
    end

    if scan_type == "delay_scan"
        filename = sample *"_"* "$surface_density_value"*"mNm-1"*"_"* polarisation_comb *"_"*split(scan)[1]*"_"*split(scan)[2] *"_"* sample_prep * dashboard_date*".h5"
    elseif scan_type == "wavenumber_scan"
        filename = sample *"_"* "$surface_density_value"*"mNm-1"*"_"* polarisation_comb *"_"*split(scan)[1]*"_"* sample_prep * dashboard_date*".h5"
    end

    
    # calculate pump wavenumber (Ekspla) for delay scan 
    ekspla_wavelength = get_metadata(raw_spectra[1])["ekspla laser"]["ekspla wavelength"][1]
    ekspla_wavenumber = round(10^7 / ekspla_wavelength, digits=2)
    
    h5open(joinpath(save_path,filename), "w") do fid
        g0 = g_create(fid, sample)
        g1 = g_create(g0, surface_density)
        g2 = g_create(g1, polarisation_comb)
        g3 = g_create(g2, scan_type)
        
        if scan_type == "delay_scan"

            if first(size(sigmatrix)) !== first(size(delay_time))
                error("Dimensions of sigmatrix and delay_time must match!!\nYou got $(first(size(sigmatrix)))-elements in sigmatrix and $(first(size(delay_time)))-elements in delay_time! ")
            elseif last(size(sigmatrix)) !== first(size(probe_wavenumbers))
                error("Dimensions of sigmatrix and probe_wavenumbers must match!!\nYou got $(last(size(sigmatrix)))-elements in sigmatrix and $(first(size(probe_wavenumbers)))-elements in probe_wavenumbers!")
            elseif first(size(sig_bleaches[1])) !== first(size(delay_time))
                error("Dimensions of sig_bleaches[i] and delay_time must match!!\nYou got $(first(size(sig_bleaches[1])))-elements in sig_bleaches[i] and $(first(size(delay_time)))-elements in delay_time! ")
            else
                

                g4 = g_create(g3, pump_resonance)
                g5 = g_create(g4, dashboard_date)
            
                g6 = g_create(g5, "Data")
                g6["sig_matrix"] = sigmatrix
                g6["ref_matrix"] = refmatrix
                g6["pump_wavenumber"] = ekspla_wavenumber
                g6["wavenumber"] = probe_wavenumbers
                g6["dltime"] = delay_time
            
                for i in 1:length(mode_name)
                    g6["sig_mean_$(mode_name[i])"] = sig_bleaches[i]
                end

                for i in 1:length(mode_name)
                    g6["ref_mean_$(mode_name[i])"] = ref_bleaches[i]
                end

                g6["comment"] = get_metadata(raw_spectra[1])["comment"][1]*add_comment
                g6["folder_name"] = directory
            end
    
            
        elseif scan_type == "wavenumber_scan"

            if pump_wavenumbers === nothing 
                pump_wavenumbers = [10^7 ./ get_metadata(raw_spectra[i])["ekspla laser"]["ekspla wavelength"][1] for i in 1:length(raw_spectra)]
            end

            if first(size(sigmatrix)) !== first(size(pump_wavenumbers))
                error("Dimensions of sigmatrix and pump_wavenumbers must match!!\nYou got $(first(size(sigmatrix)))-elements in sigmatrix and $(first(size(pump_wavenumbers)))-elements in pump_wavenumbers! ")
            elseif last(size(sigmatrix)) !== first(size(probe_wavenumbers))
                error("Dimensions of sigmatrix and probe_wavenumbers must match!!\nYou got $(last(size(sigmatrix)))-elements in sigmatrix and $(first(size(probe_wavenumbers)))-elements in probe_wavenumbers!")
            elseif first(size(sig_bleaches[1])) !== first(size(pump_wavenumbers))
                error("Dimensions of sig_bleaches[i] and pump_wavenumbers must match!!\nYou got $(first(size(sig_bleaches[1])))-elements in sig_bleaches[i] and $(first(size(pump_wavenumbers)))-elements in pump_wavenumbers! ")
            else

                g4 = g_create(g3, dashboard_date)
                
                g5 = g_create(g4, "Data")
                g5["sig_matrix"] = sigmatrix
                g5["ref_matrix"] = refmatrix
                g5["pump_wavenumber"] = pump_wavenumbers # Ekspla pump wavenumber for wavenumber scan
                g5["wavenumber"] = probe_wavenumbers
                g5["dltime"] = delay_time
                
                for i in 1:length(mode_name)
                    g5["sig_mean_$(mode_name[i])"] = sig_bleaches[i]
                end

                for i in 1:length(mode_name)
                    g5["ref_mean_$(mode_name[i])"] = ref_bleaches[i]
                end

                g5["comment"] = get_metadata(raw_spectra[1])["comment"][1]*add_comment
                g5["folder_name"] = directory    
            end
        end
    end
end

function save_dl_scan( sample::AbstractString, measurement::AbstractString,polarisation_comb::AbstractString;
                    v_surface_density::AbstractString = "SAM",
                    date = Main.date, 
                    pump_resonance::AbstractString = "" ,
                    raw_spectra = Main.raw, 
                    sigmatrix = Main.sig03, 
                    refmatrix = Main.ref03, 
                    probe_wavenumbers = Main.ν, 
                    delay_time = Main.dltime_sorted,
                    mode_name = Main.mode_name, 
                    sig_bleaches = [mean(Main.sig03[:,Main.pixel[i]], dims=2)[:,1] for i in 1:length(Main.mode_name)],
                    ref_bleaches = [mean(Main.ref03[:,Main.pixel[i]], dims=2)[:,1] for i in 1:length(Main.mode_name)],
                    add_comment= "",
                    save_path="./"
        )

    # Check if date has the right type
    if typeof(date) == String
        directory = date * "/" * foldername
        dashboard_date = date
    elseif  typeof(date) !== String 
        dashboard_date = date2dashboard(date)
        directory = dashboard_date * "/" * foldername
    elseif length(date) !== 8
        error("$date has not the type  yyyymmdd")
    end

    # sample surface density

    if v_surface_density !== "SAM"
        surface_density = "$surface_density_value mNm⁻¹"
    end
        
    # File Name and Save Path
    filename = sample*"-"*measurement*".h5"
    if save_path !== "./" 
    save_path = projectdir("data/exp_pro/Spectroscopy/sample/$filename")
    end


    
    # calculate pump wavenumber (Ekspla) for delay scan 
    ekspla_wavelength = get_metadata(raw_spectra[1])["ekspla laser"]["ekspla wavelength"][1]
    if ekspla_wavelength == 3420
        pump_resonance = "d⁻pumped"
    elseif ekspla_wavelength == 3378
        pump_resonance = "r⁻pumped"
    else
        error("""Usually we only pump d⁻ and r⁻. If you want to pump something else set the kwarg "pump_resonance" """)
    end
    ekspla_wavenumber = round(10^7 / ekspla_wavelength, digits=2)
    
    h5open(save_path, "w") do fid
        g0 = g_create(fid, sample)
        g1 = g_create(g0, surface_density)
        g2 = g_create(g1, polarisation_comb)
        g3 = g_create(g2, scan_type)
        


            if first(size(sigmatrix)) !== first(size(delay_time))
                error("Dimensions of sigmatrix and delay_time must match!!\nYou got $(first(size(sigmatrix)))-elements in sigmatrix and $(first(size(delay_time)))-elements in delay_time! ")
            elseif last(size(sigmatrix)) !== first(size(probe_wavenumbers))
                error("Dimensions of sigmatrix and probe_wavenumbers must match!!\nYou got $(last(size(sigmatrix)))-elements in sigmatrix and $(first(size(probe_wavenumbers)))-elements in probe_wavenumbers!")
            elseif first(size(sig_bleaches[1])) !== first(size(delay_time))
                error("Dimensions of sig_bleaches[i] and delay_time must match!!\nYou got $(first(size(sig_bleaches[1])))-elements in sig_bleaches[i] and $(first(size(delay_time)))-elements in delay_time! ")
            else
                

                g4 = g_create(g3, pump_resonance)
                g5 = g_create(g4, dashboard_date)
            
                g6 = g_create(g5, "Data")
                g6["sig_matrix"] = sigmatrix
                g6["ref_matrix"] = refmatrix
                g6["pump_wavenumber"] = ekspla_wavenumber
                g6["wavenumber"] = probe_wavenumbers
                g6["dltime"] = delay_time
            
                for i in 1:eachindex(mode_name)
                    g6["sig_mean_$(mode_name[i])"] = sig_bleaches[i]
                end

                for i in 1:eachindex(mode_name)
                    g6["ref_mean_$(mode_name[i])"] = ref_bleaches[i]
                end

                g6["comment"] = get_metadata(raw_spectra[1])["comment"][1]*add_comment
                g6["folder_name"] = directory
            end
    end
end
