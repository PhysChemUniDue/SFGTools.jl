const r_ssp = r"ssp"
const r_ppp = r"ppp"
const r_sps = r"sps"
const r_pss = r"pss"


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
               push!(pol_comb,match.match)
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


function get_pol_comb(s::SFSpectrum)

    comment = get_metadata(s)["comment"][1]

    get_pol_comb(comment)  
end
