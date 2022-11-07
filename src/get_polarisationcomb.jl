const r_ssp = r"ssp"
const r_ppp = r"ppp"
const r_sps = r"sps"
const r_pss = r"pss"

"""
Fetch the polarisation combination of a SFSpectrum or a direclty of a comment.\n 
\t get_pol_comb(comment::AbstractString)\n
\t get_pol_comb(s::SFSpectrum)\n


julia> get_pol_comb(raw[1])
\t "ssp"
"""
function get_pol_comb(comment::AbstractString)

    m_ssp = match(r_ssp,comment)
    m_ppp = match(r_ppp,comment)
    m_sps = match(r_sps,comment)
    m_pss = match(r_pss,comment)



       matches = [m_ssp;m_ppp;m_sps;m_pss]
    n=0
    pol_comb = []
       for match in matches
           if match !== nothing
               n += 1
               string = "$(match.match[1:3])"
               push!(pol_comb,string)
           end
       end
           if n==1
               return pol_comb[1]            
           elseif n==0
               error("The comment of the given SFSpectrum does not contain a polarization combination.")
           elseif n>1
                error("The comment of the given SFSpectrum contains the following polarisations combinations:\t
               $(pol_comb)
               ")
        
           end
end

function get_comment(s::SFSpectrum)
    comment = get_metadata(s)["comment"]

    if typeof(comment) <: Vector
        comment = comment[1]
    end
    return comment
end

function get_timestamp(s::SFSpectrum)
    timestamp = get_metadata(s)["timestamp"]

    if typeof(timestamp) <: Vector
        timestamp = timestamp[1]
    end
    return timestamp
end

function get_pol_comb(s::SFSpectrum)

    comment = get_comment(s)

    get_pol_comb(comment)  
end
