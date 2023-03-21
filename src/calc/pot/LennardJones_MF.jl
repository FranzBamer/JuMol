## lj force field
# cf. M. Falk, J.S. Langer (1998) PRE 57:7192
################################################
################################################
################################################


import Core.Array

mutable struct LJ
    Rc::Float64
    LJ(Rc) = new(Rc)
    sigma::Array{Float64}
    epsilon::Array{Float64}
    V_Rc::Vector{Float64}
    V_prime_Rc::Vector{Float64}
    optimize::Bool
end

## setting all parameters: this must be executed before the potential is used
function set_LJ_params(lj::LJ,optimize=false)
    ##
    lj.sigma = ([2*sin(pi/10)   1.0   2*sin(pi/5)])
    lj.epsilon = ([0.5   1.0   0.5])
    #### for shifting and tilding
    ## calculate potential at the cutoff
    lj.V_Rc = zeros(3)
    lj.V_Rc[1] = calc_V(lj::LJ, lj.Rc, lj.epsilon[1], lj.sigma[1]);
    lj.V_Rc[2] = calc_V(lj::LJ, lj.Rc, lj.epsilon[2], lj.sigma[2]);
    lj.V_Rc[3] = calc_V(lj::LJ, lj.Rc, lj.epsilon[3], lj.sigma[3]);
    ## calculate derivative of the potential at the cutoff
    lj.V_prime_Rc = zeros(3)
    lj.V_prime_Rc[1] = calc_V_prime(lj::LJ, lj.Rc, lj.epsilon[1], lj.sigma[1]);
    lj.V_prime_Rc[2] = calc_V_prime(lj::LJ, lj.Rc, lj.epsilon[2], lj.sigma[2]);
    lj.V_prime_Rc[3] = calc_V_prime(lj::LJ, lj.Rc, lj.epsilon[3], lj.sigma[3]);
    ## optimize speed
    lj.optimize = optimize
end

## calculate the iter-particle potential
function calc_pot(lj::LJ, type1, type2, r)
    pot = 0.0
    if r < lj.Rc
        # current id of atom1 and atom2 -> find parameters
        pos = find_pos(lj,type1,type2)
        sig = lj.sigma[pos]
        eps = lj.epsilon[pos]
        # calculate potential
        pot = calc_V(lj, r, eps, sig) - lj.V_Rc[pos] - lj.V_prime_Rc[pos]*(r-lj.Rc)
    end
  return pot
end

## calculate the inter-particle force
function calc_force(lj::LJ, type1, type2, dist_vec)
    force = zeros(3)
    if dist_vec[4] < lj.Rc
        # current id of atom1 and atom2 -> find parameters
        pos = find_pos(lj,type1,type2)
        sig = lj.sigma[pos]
        eps = lj.epsilon[pos]
        # calculate derivative
        dU_dr = calc_V_prime(lj, dist_vec[4], eps, sig)
        # find shifting correction
        du_dr_corr = lj.V_prime_Rc[pos]
        # calculate force
        force = (-1.0)*(dU_dr-du_dr_corr) * 1.0/dist_vec[4] * dist_vec[1:3]
        #correction of cutoff radius is not included
    end
    return force
end

## define pair (find the corresponding positions of the input matrices)
function find_pos(lj::LJ, type_atom1, type_atom2)
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
function calc_V(lj::LJ, r, eps, sig)
    return 4 * eps * ((sig/r)^12.0 - (sig/r)^6.0 )
    #return 4 * eps * ((sig/r)^12.0 ) #
    #return 4 * eps * (-(sig/r)^6.0 ) #
end

## calculate derivative of the potential
function calc_V_prime(lj::LJ, r, eps, sig)
    return 4 * eps * (6 * sig^6.0/r^7.0 - 12 *sig^12.0 /r^13.0)
end
