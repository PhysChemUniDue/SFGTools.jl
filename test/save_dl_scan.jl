using SFGTools
using SFGToolsScripts
function fk()
    dir = "test/sampledata/2023-01-21"
    grab(dir, getall = true)
    date = (2023,1,20)
    df = list_spectra(date = date, inexact = "DL-ODT-006", group = true)
    raw = [load_spectra(df.id[1])...]
    dark = load_spectra(df.id[2])[1]
    ν = get_ir_wavenumber(raw[1], vis = 513.35)
    mode_loc = 2875, 2937
    mode_name = ["r⁺", "rFr"]
    dltime = get_pump_delay.(raw);
    dlpos_zero = 59700;
    dltime .-= dlpos2t(dlpos_zero);
    dltime .*= -1
    perm = sortperm(dltime);
    dltime_sorted = deepcopy(dltime)
    sort!(dltime_sorted);
    raw_sorted = deepcopy(raw)
    raw_sorted = raw[perm]
    dark01 = deepcopy(dark);
    events_dark = 0
    for _= 1:10 events_dark += rm_events!(dark01,width= 10,minstd=3) end
    data01 = deepcopy(raw_sorted)
    n_events_data = zeros(Int, length(data01))
    for _= 1:10 n_events_data .+= rm_events!.(data01, width = 5, minstd = 5) end
    data02 = deepcopy(data01)
    fieldcorrection!.(data02, dark = dark01)
    data03 = sum.(data02, dims = 3)
    sig01 = [data03[i][p,1,1] for i in 1:length(data03), p in 1:size(data03[1],1)];
    ref01 = [data03[i][p,2,1] for i in 1:length(data03), p in 1:size(data03[1],1)]
    pixel = [243:272, 165:180]
    mtmp = [data03[i][j,1,1] for i in 1:length(data03), j in 1:size(data03[1],1)]
    threshold = 15
    reffac = get_reffactors(ref01, refrange = 1:10, verbose = true, threshold = threshold)
    ref02 = ref01 ./ reffac
    sig02 = sig01 ./ reffac
    ref03 = normalize(ref02, navg = 1160)
    sig03 = normalize(sig02, navg = 800)


    return date,raw,sig03,ref03,ν,dltime_sorted,mode_name,pixel
end
date,raw, sig03,ref03, ν,dltime_sorted,mode_name,pixel = fk()

save_path = "test/sampledata/exp_pro/ODT-006-006.h5"
save_dl_scan("ODT-006","006","ssp",save_path=save_path)

