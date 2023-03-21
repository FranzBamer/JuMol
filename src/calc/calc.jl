## force and potential calculations
######################################
######################################
######################################


include("pot/Yukawa_2D.jl")
include("pot/Tersoff.jl")
include("pot/LennardJones_MF.jl")
include("pot/Harm_ring_ring.jl")


mutable struct Calc
    structure
    Calc(structure) = new(structure)
    type_potential::Int64
    potential
    stress_tensor
end

## initialize type of potential before the calculation
function initialize_potential(calc::Calc, type_potential)
    calc.type_potential = type_potential
    if calc.type_potential == 7
        calc.potential = Harm_ring_ring()
        set_ring_ring_params(calc.potential)
    end
    if calc.type_potential == 13
        calc.potential = Yukawa_2D(calc.structure.rc)
        set_Yukawa_2D_params(calc.potential)
    end
    if calc.type_potential == 21
        calc.potential = Tersoff(calc.structure.rc)
        set_Tersoff_params(calc.potential,calc.structure)
    end
    if calc.type_potential == 33
        calc.potential = LJ(calc.structure.rc)
        set_LJ_params(calc.potential)
    end
    ## initialize stress tensor
    stress_tensor = zeros(3,3)
end


## set forces on all atoms to zero
function set_atom_forces_to_zero(calc::Calc)
    for i in 1:calc.structure.noa
        calc.structure.atom_list[i].force = zeros(3)
    end
    ## set also stress tensor to zero
    calc.stress_tensor = zeros(3,3)
end

## set potentials of all atoms to zero
function set_atom_potential_to_zero(calc::Calc)
    for i in 1:calc.structure.noa
        calc.structure.atom_list[i].pot = 0.0
    end
end

## calculate the sum of all pair potentials on all atoms
function calc_all_pair_pot(calc::Calc)
    set_atom_potential_to_zero(calc)
    for i in 1:calc.structure.noa
        for ii in calc.structure.atom_list[i].neighbor_indices
            pot = 0.0
            if ii > i
                if calc.type_potential == 7 || calc.type_potential == 13
                    pot = calc_pot(calc.potential,calc.structure.atom_list[i].type,
                                   calc.structure.atom_list[ii].type,
                                   calc.structure.atom_list[i].distances[4,ii-i])
                    calc.structure.atom_list[i].pot += pot
                    calc.structure.atom_list[ii].pot += pot
                end
                if calc.type_potential == 33
                    pot = calc_pot(calc.potential,calc.structure.atom_list[i].type,
                                   calc.structure.atom_list[ii].type,
                                   calc.structure.atom_list[i].distances[4,ii-i])
                    calc.structure.atom_list[i].pot += pot
                    calc.structure.atom_list[ii].pot += pot
                end

            end
            if calc.type_potential == 21
                pot = calc_pot(calc.potential,i,ii)
                calc.structure.atom_list[i].pot += pot
            end

        end
    end
end

## calculate total potential energy of the system
function calc_total_pot_en(calc::Calc)
    calc_all_pair_pot(calc)
    total_pot = 0.0
    for ii in 1:calc.structure.noa
        total_pot += calc.structure.atom_list[ii].pot
    end
    return total_pot
end

## calculate the sum of all forces on each atom
function calc_all_pair_forces(calc::Calc, calc_stress_tensor=false)
    set_atom_forces_to_zero(calc)
    for i in 1:calc.structure.noa
        for ii in calc.structure.atom_list[i].neighbor_indices
            if ii > i
                force = zeros(3)
                if calc.type_potential == 7 || calc.type_potential == 13
                    force = calc_force(calc.potential,calc.structure.atom_list[i].type,
                                       calc.structure.atom_list[ii].type,
                                       calc.structure.atom_list[i].distances[:,ii-i])
                end
                if calc.type_potential == 33
                    force = calc_force(calc.potential,calc.structure.atom_list[i].type,
                                       calc.structure.atom_list[ii].type,
                                       calc.structure.atom_list[i].distances[:,ii-i])
                end
                if calc.type_potential == 21
                    force = calc_force(calc.potential,i,ii)
                end
                calc.structure.atom_list[i].force += force
                calc.structure.atom_list[ii].force -= force
                if calc_stress_tensor == true
                    if calc.structure.atom_list[i].calc_stress == true
                        calc.stress_tensor += calc.structure.atom_list[i].distances[1:3,ii-i]*force'
                    end
                end
            end
        end
    end
    if calc_stress_tensor == true
        calc.stress_tensor = calc.stress_tensor*(-1.0)/(calc.structure.box.lx*calc.structure.box.ly*calc.structure.box.lz)
    end
end

## calculate the force for Tersoff
function calc_all_pair_forces_Tersoff(calc::Calc,calc_stress_tensor=false)
    set_atom_forces_to_zero(calc)
    for i in 1:calc.structure.noa
        calc.potential.precomputed = zeros(2,size(calc.structure.atom_list[i].neighbor_indices)[1])
        calc_b_bp_full(calc.potential,i)
        index = 1
        for ii in calc.structure.atom_list[i].neighbor_indices  # ii = j
            partial_force = zeros(3)
            b_ij = calc.potential.precomputed[1,index]
            bp_ij = calc.potential.precomputed[2,index]
            partial_force = calc_dU_i_dr_ij(calc.potential,i,ii,b_ij,bp_ij)   # dU_i_dr_ij
            calc.structure.atom_list[i].force += partial_force
            calc.structure.atom_list[ii].force -= partial_force
            index +=1
        end
    end
end


## extra function for calculating the stress tensor
function calc_stress_tensor(calc::Calc)
    calc.stress_tensor = zeros(3,3)
    #
    for i in 1:calc.structure.noa
        if calc.structure.atom_list[i].group == 0
            for ii in calc.structure.atom_list[i].neighbor_indices
                if i != ii
                    dist_vec = get_distance_between2atoms(calc.structure, i, ii)
                    force = zeros(3)
                    if calc.type_potential == 7 || calc.type_potential == 13
                        force = calc_force(calc.potential,calc.structure.atom_list[i].type,
                                           calc.structure.atom_list[ii].type,
                                           dist_vec)
                    end
                    if calc.type_potential == 33
                        force = calc_force(calc.potential,calc.structure.atom_list[i].type,
                                           calc.structure.atom_list[ii].type,
                                           dist_vec)
                    end
                    if calc.type_potential == 21
                        force = calc_force(calc.potential,i,ii)
                    end
                    #println(dist_vec[1:3])
                    #println(force')
                    #readline()
                    calc.stress_tensor += dist_vec[1:3]*force'
                end
            end
        end
    end
    calc.stress_tensor = (0.5)*calc.stress_tensor*(-1.0)/(calc.structure.box.lx*calc.structure.box.ly*calc.structure.box.lz)
end
