## Tersoff force field (Tersoff 1988, Lindsay 2010, Fan 2015)
######################################
######################################
######################################

#### https://github.com/brucefan1983/GPUMD - sample code from Zheyong Fan in C++

import Core.Array

include("../../structure/Structure.jl")


mutable struct Tersoff
    Rc::Float64
    Tersoff(Rc) = new(Rc)
    params::Array{Float64}
    structure
    precomputed
end

## setting all parameters: this must be executed before the potential is used
function set_Tersoff_params(tersoff::Tersoff,structure)
    #### parameter setting according to papers Tersoff...
    tersoff.structure = structure
    ## Lindsay (2010) Physical Review B - parameters for C-C systems
    ## A [1], B [2], λ_1 [3], λ_2 [4], λ_3 [5], n [6], c [7], β [8], d [9]  , h [10],  R [11] ,  D [12]
    ## parameter_name[column no.]
    tersoff.params = ([1393.6  346.74  3.4879  2.2119  0.0  0.72751 38049.0  1.5724e-07  4.3484  -0.57058  1.95 0.15])
    ## R+D = 2.10 A has to be given as input Rc while creating molstruc in file "tersoff_carbon.jl"
    ## λ_3 [5] = 0.0 for C-C systems, so no need to consider the exp term in ζ_ij calculation in calc_ζ
end

## calculate the inter-particle potential
function calc_pot(tersoff::Tersoff, i, j)
    pot = 0.0

    r_ij = get_distance_between2atoms(tersoff.structure,i,j)[4]
    fc_ij = calc_fc(r_ij,tersoff.params)
    fr_ij = calc_fr(r_ij,tersoff.params)
    fa_ij = calc_fa(r_ij,tersoff.params)
    b_ij  = calc_b(tersoff,i,j)

    #print("$i,$j,$b_ij")
    #print("\n")
    pot = 0.5 * fc_ij * (fr_ij - b_ij * fa_ij)
    return pot
end

## calculate the partial force term dU_i/dr_ij in force calculation [eqn. (A2) in Appendix - Fan (2015)]
function calc_dU_i_dr_ij(tersoff::Tersoff, i, j, b_ij, bp_ij)
    dU_i_dr_ij = zeros(3)
    r_ij = get_distance_between2atoms(tersoff.structure,i,j)

    # i,j terms
    fc_ij,fcp_ij = calc_fc_fcp(r_ij[4],tersoff.params)
    fr_ij,frp_ij = calc_fr_frp(r_ij[4],tersoff.params)
    fa_ij,fap_ij = calc_fa_fap(r_ij[4],tersoff.params)

    # i,j,k terms
    sum_over_k_term = zeros(3)      # terms 3 and 4 in eqn.(A2) of Fan (2015)
    sum_over_k_term = calc_ijk_terms(tersoff,i,j,fc_ij,fcp_ij,fa_ij,fap_ij, b_ij,bp_ij)
    #println(sum_over_k_term[1])
    dU_i_dr_ij = -1.0 * (0.5 * (fcp_ij * (fr_ij - b_ij * fa_ij) + fc_ij * (frp_ij - b_ij * fap_ij)) * 1.0/r_ij[4] * r_ij[1:3]
                 - 0.5 * sum_over_k_term[1:3])

    return dU_i_dr_ij
end

## calculate the smooth cutoff function value fc(r_ij) and its derivative fcp(r_ij)
function calc_fc(r,params)
    R = params[11]
    D = params[12]
    R1 = R - D
    R2 = R + D
    if r < R1
        fc = 1.0     # r < R-D
    elseif r < R2
        fc = 0.5 - 0.5 * sin(pi/2 * (r - R)/D)  #  R-D < r < R+D  (taken from Tersoff (1988))
    else
        fc = 0.0     # r > R+D
    end
    return fc
end

function calc_fc_fcp(r,params)
    R = params[11]
    D = params[12]
    R1 = R - D
    R2 = R + D
    if r < R1
        fc = 1.0     # r < R-D
        fcp = 0.0
    elseif r < R2
        fc = 0.5 - 0.5 * sin(pi/2 * (r - R)/D)  #  R-D < r < R+D
        fcp = -0.25 * pi/D * cos(pi/2 * (r - R)/D)
    else
        fc = 0.0     # r > R+D
        fcp = 0.0
    end
    return fc,fcp
end

## calculate the repulsive pair potential fr(r_ij) and its derivative frp(r_ij)
function calc_fr(r,params)
    A = params[1]
    λ1 = params[3]
    return (A * exp(-λ1*r))
end

function calc_fr_frp(r,params)
    A = params[1]
    λ1 = params[3]
    fr = A * exp(-λ1*r)
    frp = -λ1 * fr
    return fr,frp
end

## calculate the attractive pair potential fa(r_ij) and its derivative fap(r_ij)
function calc_fa(r,params)
    B = params[2]
    λ2 = params[4]
    return (B * exp(-λ2*r))
end

function calc_fa_fap(r,params)
    B = params[2]
    λ2 = params[4]
    fa = B * exp(-λ2*r)
    fap = -λ2 * fa
    return fa,fap
end

## calculate the bond angle/order term b_ij and its derivative bp_ij
function calc_b(tersoff::Tersoff,i,j)
    β = tersoff.params[8]
    n = tersoff.params[6]
    ζ_ij = calc_ζ(tersoff::Tersoff,i,j)
    b_ij = (1.0 + (β * ζ_ij)^n)^(-0.5/n)
    return b_ij
end

function calc_b_bp(tersoff::Tersoff,i,j)
    β = tersoff.params[8]
    n = tersoff.params[6]
    ζ_ij = calc_ζ(tersoff::Tersoff,i,j)
    bzn = (β * ζ_ij)^n
    bzn1 = (1.0 + bzn)
    b_ij = bzn1^(-0.5/n)
    if ζ_ij < 1e-16  ## avoid division by zero
        bp_ij = 0.0
    else
        bp_ij = -0.5 * bzn * b_ij/(bzn1 * ζ_ij)
    end
    return b_ij,bp_ij
end

function calc_b_bp_full(tersoff::Tersoff,i)
    index = 1
    for j in tersoff.structure.atom_list[i].neighbor_indices
        b_ij,bp_ij = calc_b_bp(tersoff,i,j)
        tersoff.precomputed[1,index] = b_ij
        tersoff.precomputed[2,index] = bp_ij
        index +=1
    end
end

## calculate the ζ_ij term in b_ij
function calc_ζ(tersoff::Tersoff,i,j)
    c = tersoff.params[7]
    d = tersoff.params[9]
    h = tersoff.params[10]

    ζ_ij = 0
    r_ij = get_distance_between2atoms(tersoff.structure,i,j)
    for k in tersoff.structure.atom_list[i].neighbor_indices
        if k != j
            r_ik = get_distance_between2atoms(tersoff.structure,i,k)
            fc_ik = calc_fc(r_ik[4],tersoff.params)
            cos_θ_ijk = (r_ij[1] * r_ik[1] + r_ij[2] * r_ik[2] + r_ij[3] * r_ik[3])/(r_ij[4] * r_ik[4])
            g_ijk = 1.0 + c^2/d^2 - c^2/(d^2 + (h - cos_θ_ijk)^2)
            ζ_ij += fc_ik * g_ijk
        end
    end
    return ζ_ij
end

## calculate the sum_over_k_term in force calculation (eqn. A2- Fan (2015))
function calc_ijk_terms(tersoff::Tersoff,i,j,fc_ij,fcp_ij,fa_ij,fap_ij,b_ij,bp_ij)
    c = tersoff.params[7]
    d = tersoff.params[9]
    h = tersoff.params[10]

    sum_over_k_term = zeros(3)
    index = 1
    r_ij = get_distance_between2atoms(tersoff.structure,i,j)
    d_cos_θ_ijk_dr_ij = zeros(3)
    for k in tersoff.structure.atom_list[i].neighbor_indices
        if k != j
            r_ik = get_distance_between2atoms(tersoff.structure,i,k)

            fc_ik = calc_fc(r_ik[4],tersoff.params)
            fa_ik = calc_fa(r_ik[4],tersoff.params)
            bp_ik = tersoff.precomputed[2,index]

            cos_θ_ijk = (r_ij[1] * r_ik[1] + r_ij[2] * r_ik[2] + r_ij[3] * r_ik[3])/(r_ij[4] * r_ik[4])
            g_ijk = 1.0 + c^2/d^2 - c^2/(d^2 + (h - cos_θ_ijk)^2)
            gp_ijk = -2.0 * c^2 * (h - cos_θ_ijk)/((d^2 + (h - cos_θ_ijk)^2)^2)

            d_cos_θ_ijk_dr_ij = 1.0/r_ij[4] * (1.0/r_ik[4] * r_ik[1:3] - 1.0/r_ij[4] * r_ij[1:3] * cos_θ_ijk)

            sum_over_k_term +=  fc_ik * fcp_ij * fa_ik * bp_ik * g_ijk * 1.0/r_ij[4] * r_ij[1:3] +
                                fc_ik * fc_ij * gp_ijk * d_cos_θ_ijk_dr_ij[1:3] * (fa_ij * bp_ij + fa_ik * bp_ik)
        end
        index += 1
    end
    return sum_over_k_term
end

## calculate the g_ijk term
function calc_g_ijk(tersoff::Tersoff,i,j,k)
    c = tersoff.params[7]
    d = tersoff.params[9]
    h = tersoff.params[10]

    r_ij = get_distance_between2atoms(tersoff.structure,i,j)
    r_ik = get_distance_between2atoms(tersoff.structure,i,k)

    cos_θ_ijk = (r_ij[1] * r_ik[1] + r_ij[2] * r_ik[2] + r_ij[3] * r_ik[3])/(r_ij[4] * r_ik[4])
    g_ijk = 1.0 + c^2/d^2 - c^2/(d^2 + (h - cos_θ_ijk)^2)

    return(cos_θ_ijk,g_ijk)
end
