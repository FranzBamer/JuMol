## Harmonic ring-ring-potential
######################################
######################################
######################################
#
# atom type from 4 to 10
#
#
#
# (c) Franz Bamer, Nov-2020
######################################


import Core.Array

mutable struct Harm_ring_ring
    kh::Float64
    bl::Float64
    Harm_ring_ring(kh=1.0,bl=1.0) = new(kh,bl)
    R_0_mat::Array{Float64}
end



function calc_r_0(hrr::Harm_ring_ring,type1,type2)
    term1 = 1.0/(tan(pi/float(type1)));
    term2 = 1.0/(tan(pi/float(type2)));
    return (term1+term2)*0.5;
end

## setting all parameters: this must be executed before the potential is used
function set_ring_ring_params(hrr::Harm_ring_ring)
    #### parameter setting according to Morley, 2018
    hrr.R_0_mat = zeros(10,10)
    hrr.R_0_mat[1,1] = hrr.bl # in case of a standard harmonic potential
    for i=4:10
        for ii=4:10
            hrr.R_0_mat[i,ii] = calc_r_0(hrr,i,ii)*hrr.bl
        end
    end
    #pretty_table(hrr.R_0_mat,
    #	tf = borderless, noheader = true, crop = :horizontal, formatters = ft_round(3))
end

## calculate the inter-particle potential
function calc_pot(hrr::Harm_ring_ring, type1, type2, r)
    r_0 = hrr.R_0_mat[type1,type2]
    pot = calc_V(hrr, r, r_0)
  return pot
end

## calculate the inter-particle force
function calc_force(hrr::Harm_ring_ring, type1, type2, dist_vec)
        r_0 = hrr.R_0_mat[type1,type2]
        dU_dr = calc_V_prime(hrr, dist_vec[4], r_0)
        force = (-1.0)*(dU_dr) * 1.0/dist_vec[4] * dist_vec[1:3]
        #println(force)
        #readline()
    return force
end


## calculate potential
function calc_V(hrr::Harm_ring_ring, r, r_0)
    return hrr.kh * 1.0/(r_0) * (r-r_0) * (r-r_0);
end

## calculate derivative of the potential
function calc_V_prime(hrr::Harm_ring_ring, r, r_0)
    return hrr.kh * 1.0/(r_0) * 2.0*(r-r_0)
end
