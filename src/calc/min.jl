## minimizer
######################################
######################################
######################################

include("int.jl")

## struct Min definition
mutable struct Min
    structure
    calc
    alpha_min::Float64
    full_update_num::Int64
    acc_factor::Float64
    output_filename::String
    Min(structure,calc,
        alpha_min=1.0e-2,
        full_update_num=10,
        acc_factor=1.01,
        output_filename="None") = new(structure,calc,
                                  alpha_min,
                                  full_update_num,
                                  acc_factor,
                                  output_filename)
    num_steps_max::Int64
    num_steps_real::Int64
    step_list::Vector{Int64}
    pot_en_hist::Vector{Float64}
    force_hist::Vector{Float64}
end

## run steepest decent algorithm
function run_sd(min::Min, tolerance, num_steps, group_list=[0], full_update=true, save_lammpstrj=false, save_step=100)
    ##
    if save_lammpstrj == true
        open_file(min.structure.reader,"min_hist_output.lammpstrj", "w")
    end
    min.num_steps_max = num_steps
    min.step_list = zeros(num_steps)
    min.pot_en_hist = zeros(num_steps)
    min.force_hist = zeros(num_steps)
    ##
    max_disp_step = 0.0
    ## full update before the minimization starts
    if full_update==true
        if min.structure.linked_cells_bool == true
            update_cell(structure)
        	linked_cell_list(structure)
        	construct_neigh_cell_list(structure)
        end
        update_distances(min.structure)
        update_neighbor_lists(min.structure)
    end
    ## control vector for tolerance
    contr_vec = zeros(min.structure.noa*3)
    cntr = 1
    ## minimization loop
    for i in 1:min.structure.noa
        for ii in length(group_list)
            if min.structure.atom_list[i].group == group_list[ii]
                contr_vec[cntr] = 1
                contr_vec[cntr+1] = 1
                contr_vec[cntr+2] = 1
            end
        end
        cntr += 3
    end
    ## run through the number of minimization steps
    for i in 1:num_steps
        ## step counter
        min.step_list[i] = i
        ## force calculation
        calc_all_pair_forces(min.calc)
        ## add the step downwards to the position of every atom
        max_delta_step = 0.0
        for ii in 1:min.structure.noa
            for iii in length(group_list)
                if min.structure.atom_list[ii].group == group_list[iii]
                    delta_step = min.structure.atom_list[ii].force*min.alpha_min
                    min.structure.atom_list[ii].pos += delta_step
                    len_step = norm(delta_step)
                    if len_step > max_delta_step
                        max_delta_step = len_step
                    end
                end
            end
        end
        ## update distances
        max_disp_step += max_delta_step
        if max_disp_step > min.structure.rskin*0.5
            if full_update==true
                if min.structure.linked_cells_bool == true
                    update_cell(structure)
                	linked_cell_list(structure)
                	construct_neigh_cell_list(structure)
                end
                update_distances(min.structure)
                update_neighbor_lists(min.structure)
                max_disp_step = 0.0
            else
                update_distances_partly(min.structure)
            end
        else
            update_distances_partly(min.structure)
        end
        ## potential calculation
        U = calc_total_pot_en(min.calc)
        if save_lammpstrj == true
            if (i-1)%save_step == 0
                write_box(min.structure.reader,min.structure,i)
            end
        end
        #save_pot_en(min, i, U)
        ## force_tolerance achieved -> break
        F = get_global_force_vec(min.structure) .* contr_vec
        f_scal = dot(F,F)
        #save_force(min, i, f_scal)
        print("\e[2K")
        print("\r potential: ", U, " force norm: ", f_scal, " alpha_min: ", min.alpha_min)
        if f_scal < tolerance
            println("\nforce tolerance achieved at step: ", i)
            min.num_steps_real = i
            break
        end
        ## update acutal number of steps done
        min.num_steps_real = i
    end
    if save_lammpstrj == true
        close_file(min.structure.reader)
    end
end

## get the potential energy and save it
function save_pot_en(min::Min, step_num, pot_en)
    if step_num == 1
        min.pot_en_hist = zeros(min.num_steps_max)
    end
    min.pot_en_hist[step_num] = pot_en
end

## get the force and save it
function save_force(min::Min, step_num, f_scal)
    min.force_hist[step_num] = f_scal
end

## run conjugate gradient algorithm
function run_cg(min::Min,num_steps,tolerance = 1.0e-4, group_list=[0],full_update=true, save_lammpstrj=false, save_step=100)
    ##
    if save_lammpstrj == true
        open_file(min.structure.reader,"min_hist_output_cg.lammpstrj", "w")
    end
    ##
    min.num_steps_max = num_steps
    min.step_list = zeros(num_steps)
    min.pot_en_hist = zeros(num_steps)
    min.force_hist = zeros(num_steps)
    ## full update before the minimization starts
    max_disp_step = 0.0
    if full_update
        if min.structure.linked_cells_bool == true
            update_cell(min.structure)
            linked_cell_list(min.structure)
            construct_neigh_cell_list(min.structure)
        end
        update_distances(min.structure)
        update_neighbor_lists(min.structure)
    end
    ## minimization contribution vector (minimize only those atoms that are part of a certain group)
    println("number of atoms:")
    println(min.structure.noa)
    contr_vec = zeros(min.structure.noa*3)
    cntr = 1
    for i in 1:min.structure.noa
        for ii in length(group_list)
            if min.structure.atom_list[i].group == group_list[ii]
                contr_vec[cntr] = 1
                contr_vec[cntr+1] = 1
                contr_vec[cntr+2] = 1
            end
        end
        cntr += 3
    end
    c1 = 0.0001
    c2 = 0.1
    ## initial calculation
    U = calc_total_pot_en(min.calc)
    Up1 = 0.0
    r = get_global_atom_pos_vec(min.structure)
    rp1 = zeros(min.structure.noa*3)
    calc_all_pair_forces(min.calc)
    F = get_global_force_vec(min.structure) .* contr_vec
    Fp1 = zeros(min.structure.noa*3) .* contr_vec
    p = get_global_force_vec(min.structure) .* contr_vec
    ## run minimization loop
    for i in 1:num_steps
        ## step counter
        min.step_list[i] = i
        ## testing alpha using the WOLFE conditions
        WOLFE_TOTAL = false
        while WOLFE_TOTAL == false
            WOLFE1 = false
            WOLFE2 = false
            delta_step_vec = min.alpha_min*p
            rp1 = r + delta_step_vec
            ## update positions, calculate potential engergy and forces
            set_global_pos_vec(min.structure,rp1)
            #if i%10 == 0
            if max_disp_step > min.structure.rskin*0.5 # half of the skin distance because center atom can move in oposite direction
                if full_update==true
                    if min.structure.linked_cells_bool == true
                        update_cell(min.structure)
                    	linked_cell_list(min.structure)
                    	construct_neigh_cell_list(min.structure)
                    end
                    update_distances(min.structure)
                    update_neighbor_lists(min.structure)
                    max_disp_step = 0.0
                else
                    update_distances_partly(min.structure)
                end
            else
                update_distances_partly(min.structure)
            end
            Up1 = calc_total_pot_en(min.calc)
            calc_all_pair_forces(min.calc)
            Fp1 = get_global_force_vec(min.structure) .* contr_vec
            if Up1 <= U + c1*min.alpha_min*dot(p,F*(-1.0))
                WOLFE1 = true
            else
                WOLFE1 = false
            end
            if dot(Fp1*(-1.0),p) <= c2*dot(F*(-1.0),p)
                WOLFE2 = true
            else
                WOLFE2 = false
            end
            ## check possibilities and adapt alpha_min accordingly
            if WOLFE1 == true && WOLFE2 == true
                WOLFE_TOTAL = true
                ## Polak-Ribiere for the next step
                beta = dot(Fp1,Fp1-F) / dot(F,F)
                p = Fp1 + beta*p
                r = rp1
                U = Up1
                F = Fp1
                min.alpha_min = min.alpha_min*min.acc_factor
                ## check maximum distance of the step and add it to the value
                max_disp = 0.0
                for ii in 1:min.structure.noa
                    max_disp_check = ( delta_step_vec[3*ii-2]*delta_step_vec[3*ii-2] +
                                       delta_step_vec[3*ii-1]*delta_step_vec[3*ii-1] +
                                       delta_step_vec[3*ii-0]*delta_step_vec[3*ii-0]   )
                    if max_disp < max_disp_check
                        max_disp = max_disp_check
                    end
                end
                max_disp_step += sqrt(max_disp)
                ## save the potential energy
                save_pot_en(min, i, U)
                ## save the scalar force
                save_force(min, i, dot(F,F))
                ## save lammpstrj if necessary
                if save_lammpstrj == true
                    if (i-1)%save_step == 0
                        write_box(min.structure.reader,min.structure,i)
                    end
                end
            else
                WOLFE_TOTAL = false
                #println("step: ", i)
                #println("WOLFE1: ", WOLFE1)
                #println("WOLFE2: ", WOLFE2)
                ## decrease step size
                min.alpha_min = min.alpha_min*0.5
                #println("alpha_min: ", min.alpha_min)
                ## search alpha too small -> break
                if min.alpha_min < 1.0e-4
                    println("break criterion alpha_min (inside): ", i)
                    ## update acutal number of steps done
                    min.num_steps_real = i
                    break
                end
            end
        end
        ## search alpha too small -> break
        if min.alpha_min < 1.0e-4
            println("break criterion alpha_min: ", i)
            ## update acutal number of steps done
            min.num_steps_real = i
            break
        end
        ## force_tolerance achieved -> break
        print("\e[2K")
        print("\r potential: ", U, " force norm: ", dot(F,F), " alpha_min: ", min.alpha_min)
        if dot(F,F) < tolerance
            println("\nforce tolerance achieved at step: ", i)
            ## update acutal number of steps done
            min.num_steps_real = i
            break
        end
        ## update acutal number of steps done
        min.num_steps_real = i
    end
    if save_lammpstrj == true
        close_file(min.structure.reader)
    end
end




function nvt_steps(min::Min)
    if min.structure.linked_cells_bool == true
        update_cell(min.structure)
        linked_cell_list(min.structure)
        construct_neigh_cell_list(min.structure)
    end
    update_distances(min.structure)
    update_neighbor_lists(min.structure)
    #
    num_steps = 100
    molunits = Jumol.Units(0)
    Int = Integ(min.structure, min.calc, molunits)
    run_nvt(Int, num_steps, 0.5, 0.5, false, 2, false)
end
