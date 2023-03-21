## calculation of Hessian and related quantities
using LinearAlgebra

mutable struct Hessian
    structure
    calc
    Hessian(structure,calc) = new(structure,calc)
    hessian_mat
    evals
    evecs
    ##
    hessian_mat_sparse
end

## main function to calculate Hessian
function calc_hessian(hessian::Hessian, delta)
    hessian.hessian_mat = zeros(2*hessian.structure.noa,2*hessian.structure.noa)

    total_sub_matrices = Int(0.5 * hessian.structure.noa * (hessian.structure.noa - 1) +
                         hessian.structure.noa)

    cntr = 0

    for row in 1:hessian.structure.noa  # row -> pos_atom
        for col in 1:row                # col -> force_atom
            cntr += 1
            # println("Calculating hessian : $cntr / $total_sub_matrices")
            if hessian.calc.type_potential == 33 #FIXME potential must be changed back to 21
                calc_hessian_block_matrix(hessian, row, col, delta)
            elseif row == col || row in hessian.structure.atom_list[col].neighbor_indices
                calc_hessian_block_matrix(hessian, row, col, delta)
            else
                # println("skip")
            end
        end
    end
    ## symmetric entries
    for row in 1:hessian.structure.noa*2
        for col in (row+1):hessian.structure.noa*2
            hessian.hessian_mat[row,col] = hessian.hessian_mat[col,row]
        end
    end
    ##
    hessian.hessian_mat = hessian.hessian_mat * (-1.0)
end

## calculate hessian block matrices
function calc_hessian_block_matrix(hessian::Hessian, row, col, delta)
    # F_x, pos direction-x (11),  F_y, pos direction-x(21)
    force_der_11,force_der_21 = calc_force_derivative_hessian(hessian,row,1,col,delta)
    hessian.hessian_mat[2*row-1,2*col-1] = force_der_11  # left top
    hessian.hessian_mat[2*row,2*col-1] =  force_der_21   # left bottom
    # F_x, pos direction-y(12) ,  F_y, pos direction-y(22)
    force_der_12, force_der_22 = calc_force_derivative_hessian(hessian,row,2,col,delta)
    hessian.hessian_mat[2*row-1,2*col] =  force_der_12   # right top
    hessian.hessian_mat[2*row,2*col] = force_der_22      # right bottom
end

## function for calculating force derivative
function calc_force_derivative_hessian(hessian::Hessian,pos_atom_index,pos_dir,force_atom_index,delta)
    ## move atom in positive direction by delta
    hessian.structure.atom_list[pos_atom_index].pos[pos_dir] += delta
    update_distances_hessian(hessian, force_atom_index, pos_atom_index)
    #
    hessian.calc.structure.atom_list[force_atom_index].force = zeros(3) # set_atom_forces_to_zero for "force_atom"
    calc_atom_force_hessian(hessian,force_atom_index, pos_atom_index)
    #
    f_dir_plus = hessian.calc.structure.atom_list[force_atom_index].force #[force_dir]

    ## move atom in negative direction by delta
    hessian.structure.atom_list[pos_atom_index].pos[pos_dir] -= 2.0 * delta
    update_distances_hessian(hessian, force_atom_index, pos_atom_index)
    #
    hessian.calc.structure.atom_list[force_atom_index].force = zeros(3) # set_atom_forces_to_zero for "force_atom"
    calc_atom_force_hessian(hessian,force_atom_index, pos_atom_index)
    #
    f_dir_minus = hessian.calc.structure.atom_list[force_atom_index].force #[force_dir]

    ## move atom back to its initial position and re-calcuate initial force (why?)
    hessian.structure.atom_list[pos_atom_index].pos[pos_dir] += delta
    ## may_be_update distances
    update_distances_hessian(hessian, force_atom_index, pos_atom_index)

    f_derivative = (f_dir_plus[1:2]-f_dir_minus[1:2])/(2.0 * delta)

    return f_derivative[1],f_derivative[2]
end

## function for calculating force on atom with index i
function calc_atom_force_hessian(hessian::Hessian,i, j)
    if hessian.calc.type_potential == 21
        cal_three_body_atom_force_hessian(hessian, i)
    else
        # cal_two_body_atom_force_hessian(hessian, i)
        cal_two_body_atom_force_hessian_optimized(hessian::Hessian, i, j)
    end
end

## function for calculating two-body force on atom with index 'i'
function cal_two_body_atom_force_hessian(hessian::Hessian,i)
    for ii in hessian.structure.atom_list[i].neighbor_indices
        force = zeros(3)
        if hessian.calc.type_potential == 13
            force = calc_force(hessian.calc.potential,hessian.calc.structure.atom_list[i].type,
                               hessian.calc.structure.atom_list[ii].type,
                               get_distance_between2atoms(hessian.structure,i,ii))
            hessian.calc.structure.atom_list[i].force += force
        end
        if hessian.calc.type_potential == 33
            force = calc_force(hessian.calc.potential,hessian.calc.structure.atom_list[i].type,
                               hessian.calc.structure.atom_list[ii].type,
                               get_distance_between2atoms(hessian.structure,i,ii))
            hessian.calc.structure.atom_list[i].force += force
        end
    end
end

## calculate two body atom force optimized
function cal_two_body_atom_force_hessian_optimized(hessian::Hessian,i,j)
    if j == i
        for ii in hessian.structure.atom_list[i].neighbor_indices
            force = zeros(3)
            if hessian.calc.type_potential == 13
                force = calc_force(hessian.calc.potential,hessian.calc.structure.atom_list[i].type,
                                   hessian.calc.structure.atom_list[ii].type,
                                   get_distance_between2atoms(hessian.structure,i,ii))
                hessian.calc.structure.atom_list[i].force += force
            end
            if hessian.calc.type_potential == 33
                force = calc_force(hessian.calc.potential,hessian.calc.structure.atom_list[i].type,
                                   hessian.calc.structure.atom_list[ii].type,
                                   get_distance_between2atoms(hessian.structure,i,ii))
                hessian.calc.structure.atom_list[i].force += force
            end
        end
    else
        force = zeros(3)
        if hessian.calc.type_potential == 13
            force = calc_force(hessian.calc.potential,hessian.calc.structure.atom_list[i].type,
                               hessian.calc.structure.atom_list[j].type,
                               get_distance_between2atoms(hessian.structure,i,j))
            hessian.calc.structure.atom_list[i].force += force
        end
        if hessian.calc.type_potential == 33
            force = calc_force(hessian.calc.potential,hessian.calc.structure.atom_list[i].type,
                               hessian.calc.structure.atom_list[j].type,
                               get_distance_between2atoms(hessian.structure,i,j))
            hessian.calc.structure.atom_list[i].force += force
        end
    end
end

## function for calculating three-body force (Tersoff) on atom with index i
function cal_three_body_atom_force_hessian(hessian::Hessian,i)
    if hessian.calc.type_potential == 21
        ## calculate the sum of terms -> Σ dU_i_dr_ij
        hessian.calc.potential.precomputed = zeros(2,size(hessian.calc.structure.atom_list[i].neighbor_indices)[1])
        calc_b_bp_full(hessian.calc.potential,i)
        counter = 1
        for ii in hessian.calc.structure.atom_list[i].neighbor_indices  # ii = j
            partial_force = zeros(3)
            b_ij = hessian.calc.potential.precomputed[1,counter]
            bp_ij = hessian.calc.potential.precomputed[2,counter]
            partial_force = calc_dU_i_dr_ij(hessian.calc.potential,i,ii,b_ij,bp_ij)   # dU_i_dr_ij term
            hessian.calc.structure.atom_list[i].force += partial_force
            counter +=1
        end
        ## calculate the sum of terms -> -Σ dU_j_dr_ji
        for ii in hessian.calc.structure.atom_list[i].neighbor_indices
            hessian.calc.potential.precomputed = zeros(2,
                                                size(hessian.calc.structure.atom_list[ii].neighbor_indices)[1])
            calc_b_bp_full(hessian.calc.potential,ii)
            counter = findfirst(isequal(i), hessian.calc.structure.atom_list[ii].neighbor_indices)
            partial_force = zeros(3)
            b_ij = hessian.calc.potential.precomputed[1,counter]
            bp_ij = hessian.calc.potential.precomputed[2,counter]
            partial_force = calc_dU_i_dr_ij(hessian.calc.potential,ii,i,b_ij,bp_ij)   # dU_j_dr_ji term
            hessian.calc.structure.atom_list[i].force -= partial_force
        end
    end
end

#### eigenvalue eigenvector calcuation
## function for calculating eigen values and eigen vectors
function solve_eig(hessian::Hessian, calc_eigen_vector = false, num_eigen_vectors = -1)
    if num_eigen_vectors == -1
        ef = eigen(Symmetric(hessian.hessian_mat))  # all eigen-values/eigen-vectors
    elseif num_eigen_vectors > 0
        ef = eigen(Symmetric(hessian.hessian_mat), 1:num_eigen_vectors)   #k smallest eigenvalues/vectors
    else
        print("Error while calling solve_eig")
    end
    hessian.evals = ef.values
    if calc_eigen_vector == true
        hessian.evecs = ef.vectors
    end
end

## update distance functions specific to two body force Hessian calculation
function update_distances_hessian(hessian::Hessian, force_atom_index, pos_atom_index)
    if hessian.calc.type_potential == 21
        # update_distances(hessian.structure)
        update_distances_partly(hessian.structure)
    else
        # update_distances(hessian.structure)
        updates_distances_hessian_two_body_force(hessian.structure, force_atom_index, pos_atom_index)
    end
end

function updates_distances_hessian_two_body_force(structure::Structure, i, j)
    if j == i
        for ii in structure.atom_list[i].neighbor_indices
             update_distances_between_two_atoms(structure::Structure, i, ii)
         end
    else
        update_distances_between_two_atoms(structure::Structure, i, j)
    end
end

function update_distances_between_two_atoms(structure::Structure, i, ii)
	if ii > i
		dist_vec = calc_distance_betw2atoms(structure, i, ii)
		structure.atom_list[i].distances[:,ii-i] = dist_vec
    else
        dist_vec = calc_distance_betw2atoms(structure, ii, i)
		structure.atom_list[ii].distances[:,i-ii] = dist_vec
	end
end
