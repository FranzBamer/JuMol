## Yukawa force field (Roy, 2018)
######################################
######################################
######################################


import Core.Array

mutable struct Yukawa_2D
    Rc::Float64
    Yukawa_2D(Rc) = new(Rc)
    params::Array{Float64}
    kappa::Float64
    V_Rc::Vector{Float64}
    V_prime_Rc::Vector{Float64}
    optimize::Bool
end

## setting all parameters: this must be executed before the potential is used
function set_Yukawa_2D_params(yukawa::Yukawa_2D,optimize=false)
    #### parameter setting according to Roy, 2018
    ## row 1: sigma
    ## row 2: q
    ## col 1: Si-Si , col 2: Si-O, col 3: O-O
    yukawa.params = ([2.250  1.075  0.900;
                      1.500 -1.000  0.670])
    #### long range parameter kappa
    ## Carre et al. (2007) The journal of chemical physics 127:114512.
    yukawa.kappa = 1.0/5.649
    #### for shifting and tilding
    ## calculate potential at the cutoff
    yukawa.V_Rc = zeros(3)
    yukawa.V_Rc[1] = calc_V(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1,1],yukawa.params[2,1]); #sigma(Si-Si), q(Si-Si)
    yukawa.V_Rc[2] = calc_V(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1,2],yukawa.params[2,2]); #Sigma(Si-O) , q(Si-O)
    yukawa.V_Rc[3] = calc_V(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1,3],yukawa.params[2,3]); #Sigma(O-O)  , q(O-O)
    ## calculate derivative of the potential at the cutoff
    yukawa.V_prime_Rc = zeros(3)
    yukawa.V_prime_Rc[1] = calc_V_prime(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1,1],yukawa.params[2,1]); #sigma(Si-Si), q(Si-Si)
    yukawa.V_prime_Rc[2] = calc_V_prime(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1,2],yukawa.params[2,2]); #Sigma(Si-O) , q(Si-O)
    yukawa.V_prime_Rc[3] = calc_V_prime(yukawa::Yukawa_2D,yukawa.Rc,yukawa.params[1,3],yukawa.params[2,3]); #Sigma(O-O)  , q(O-O)
    ## optimize speed
    yukawa.optimize = optimize
end

## calculate the inter-particle potential
function calc_pot(yukawa::Yukawa_2D, type1, type2, r)
    pot = 0.0
    if r < yukawa.Rc
        # current id of atom1 and atom2 -> find parameters
        pos = find_pos(yukawa,type1,type2)
        sig = yukawa.params[1,pos]
        q = yukawa.params[2,pos]
        # calculate potential
        if yukawa.optimize
            pot = calc_V(yukawa,r, sig, q) - yukawa.V_Rc[pos] - yukawa.V_prime_Rc[pos]*(r-yukawa.Rc)
        else
            pot = calc_V_opt(yukawa,r, sig, q) - yukawa.V_Rc[pos] - yukawa.V_prime_Rc[pos]*(r-yukawa.Rc)
        end
    end
  return pot
end

## calculate the inter-particle force
function calc_force(yukawa::Yukawa_2D, type1, type2, dist_vec)
    force = zeros(3)
    if dist_vec[4] < yukawa.Rc
        # current id of atom1 and atom2 -> find parameters
        pos = find_pos(yukawa,type1,type2)
        sig = yukawa.params[1,pos]
        q = yukawa.params[2,pos]
        # calculate derivative
        if yukawa.optimize
            ## calculate the inverse of r
            one_over_r = 1.0/dist_vec[4]
            dU_dr = calc_V_prime_opt(yukawa, dist_vec[4], one_over_r, sig, q)
        else
            dU_dr = calc_V_prime(yukawa, dist_vec[4], sig, q)
        end
        # find shifting correction
        du_dr_corr = yukawa.V_prime_Rc[pos]
        # calculate force
        force = (-1.0)*(dU_dr - du_dr_corr) * 1.0/dist_vec[4] * dist_vec[1:3]
    end
    return force
end

## define pair (find the corresponding positions of the input matrices)
function find_pos(yukawa::Yukawa_2D, type_atom1, type_atom2)
    col = 1
    # Si-Si
    if type_atom1 == 1 && type_atom2 == 1
        col = 1
    end
    # Si-O
    if type_atom1 == 1 && type_atom2 == 2
        col = 2
    end
    # O-Si
    if type_atom1 == 2 && type_atom2 == 1
        col = 2
    end
    # O-O
    if type_atom1 == 2 && type_atom2 == 2
        col = 3
    end
    return col
end

## calculate potential
function calc_V(yukawa::Yukawa_2D, r, sig, q)
    return (sig/r)^12.0 + q/r * exp(-yukawa.kappa*r)
end

## calculate derivative of the potential
function calc_V_prime(yukawa::Yukawa_2D, r, sig, q)
    return (-12.0)/sig*(sig/r)^13.0 - yukawa.kappa*q/r * exp(-yukawa.kappa*r) - q/(r*r) * exp(-yukawa.kappa*r)
end

## calculate potential (optimized version)
function calc_V_opt(yukawa::Yukawa_2D, r, sig, q)
    one_over_r = 1.0/r
    sig_over_r = one_over_r*sig
    sig_over_r2 = sig_over_r*sig_over_r
    sig_over_r4 = sig_over_r2*sig_over_r2
    sig_over_r8 = sig_over_r4*sig_over_r4
    q_over_r = one_over_r*q
    return sig_over_r8*sig_over_r4 + q_over_r*exp(-yukawa.kappa*r)
end

## calculate derivative of the potential (optimized version)
function calc_V_prime_opt(yukawa::Yukawa_2D, r, one_over_r, sig, q)
    q_over_r = one_over_r*q
    one_over_r2 = one_over_r*one_over_r
    sig_over_r2 = one_over_r2*sig*sig
    sig_over_r4 = sig_over_r2*sig_over_r2
    sig_over_r8 = sig_over_r4*sig_over_r4
    q_over_r2 = one_over_r2*q
    return (-12.0)*sig_over_r8*sig_over_r4*one_over_r - yukawa.kappa*q_over_r * exp(-yukawa.kappa*r) - q_over_r2 * exp(-yukawa.kappa*r)
end
