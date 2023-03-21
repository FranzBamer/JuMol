## time integrators
######################################
######################################
######################################

using Random, Distributions

mutable struct Integ
    structure
    calc
	units
    delta_t::Float64
	output_filename::String
	steps_save::Int64
    Integ(structure,calc,
	      units, delta_t=1.0e-2,
		  output_filename="None",
		  steps_save=1) = new(structure,calc,
		                                units,delta_t,
										output_filename,
										steps_save)
	## time step vector
	timestep_vec
	## NVT vector
	temp_vec_target
	temp_vec
	zeta_vec
end


## NVE ensemble
function run_nve(int::Integ, num_steps)
	## write file if necessary
	if int.output_filename != "None"
		open_file(int.structure.reader, int.output_filename, "w")
	end
	## starting NVE integration
	println("####")
	println("starting NVE integration...")
	## calculate initial force
	calc_all_pair_forces(int.calc)
	## get storage for the calculation parameters
	# velocity at the half step
	vel_half = zeros(int.structure.noa,3)
	# maximum walking distance
	max_dist = 0
	# maximum vel at each time step
	max_vel = 0
 	for i in 1:num_steps
		print("\e[2K")
 		print("\rstep: ", i)
		# update position and calculate intermediate velocity
        for ii in 1:int.structure.noa
            vel_half[ii,:] = int.structure.atom_list[ii].vel + 0.5 * int.structure.atom_list[ii].force/int.structure.atom_list[ii].mass * int.units.force_factor*int.delta_t
            int.structure.atom_list[ii].pos += vel_half[ii,:] * int.delta_t
        end
        # update distances
		for ii in 1:int.structure.noa
			v = norm(int.structure.atom_list[ii].vel)
			if v > max_vel
				max_vel = v
			end
		end
		max_dist += max_vel * int.delta_t
		max_vel = 0

		if max_dist > int.structure.rskin * 0.5
			if int.structure.linked_cells_bool == true
				update_cell(int.structure)
				linked_cell_list(int.structure)
				construct_neigh_cell_list(int.structure)
			end
	    	update_distances(int.structure)
			update_neighbor_lists(int.structure)
			max_dist = 0
		else
	    	update_distances_partly(int.structure)
		end
        # update force
        calc_all_pair_forces(int.calc)
        # update velocity
        for ii in 1:int.structure.noa
	    	int.structure.atom_list[ii].vel = vel_half[ii,:] + 0.5 * int.structure.atom_list[ii].force/int.structure.atom_list[ii].mass * int.units.force_factor* int.delta_t
		end
		put_atoms_back_to_box(int.structure)
		if int.output_filename != "None"
			if i%int.steps_save == 0
				Jumol.write_box(int.structure.reader,int.structure,i)
			end
		end
    end
	## close output file
	if int.output_filename != "None"
		close_file(int.structure.reader)
	end
end

## NVT ensemble
function run_nvt(int::Integ, num_steps, t_beg, t_end, tau, order_of_problem=2, zeta=0.0)
	## write file if necessary
	if int.output_filename != "None"
		open_file(int.structure.reader, int.output_filename, "w")
	end
	## set atom velocities
	set_atom_velocities(int, 2, t_beg)
	## get temperature history
	int.temp_vec_target= get_temp_vec(int,t_beg,t_end,num_steps)
	int.timestep_vec = LinRange(1,num_steps,num_steps)
	## starting NVT integration
	println("####")
	println("starting NVT integration...")
	## calculate initial force
	calc_all_pair_forces(int.calc)
	## get storage for the calculation parameters
	# velocity at the half step
	vel_half = zeros(int.structure.noa,3)
	# maximum walking distance
	max_dist = 0
	# maximum vel at each time step
	max_vel = 0
	# target temperature (must be defined as input parameter)
	T_tar = 0.0
	# initialize actual temperature
	T = 0.0
	# initialize kinetic energy
	KE = 0
	# degrees of freedom for 2d
	dof = order_of_problem * int.structure.noa
	# temperature vector at all time steps
	int.temp_vec = zeros(num_steps)
	int.zeta_vec = zeros(num_steps)
	for i in 1:num_steps
		## temperature
		T_tar = int.temp_vec_target[i]
		print("\e[2K")
        print("\rstep: ", i)
		#update position and calculate intermediate velocity
    	for ii in 1:int.structure.noa
        	vel_half[ii,:] = int.structure.atom_list[ii].vel + 0.5 * (int.structure.atom_list[ii].force /int.structure.atom_list[ii].mass* int.units.force_factor - zeta * int.structure.atom_list[ii].vel)* int.delta_t
			int.structure.atom_list[ii].pos += vel_half[ii,:] * int.delta_t
		end
		# update distances
		for ii in 1:int.structure.noa
			v = norm(int.structure.atom_list[ii].vel)
			if v > max_vel
				max_vel = v
			end
		end
		max_dist += max_vel * int.delta_t
		max_vel = 0
		if max_dist > int.structure.rskin*0.5
			if int.structure.linked_cells_bool == true
                update_cell(int.structure)
            	linked_cell_list(int.structure)
            	construct_neigh_cell_list(int.structure)
            end
    		update_distances(int.structure)
			update_neighbor_lists(int.structure)
			max_dist = 0
		else
    		update_distances_partly(int.structure)
		end
    	#update force
    	calc_all_pair_forces(int.calc)
		# calculate kinetic energy
		KE = 0
		for ii in 1:int.structure.noa
			KE += 0.5 * int.structure.atom_list[ii].mass * dot(int.structure.atom_list[ii].vel,int.structure.atom_list[ii].vel)
		end
		# calculate temperature
		T = 2 * KE / (dof * int.units.Kb )
		zeta += 0.5 * int.delta_t * (T/T_tar - 1.0)/tau^2
		KE = 0
		for ii in 1:int.structure.noa
			KE += 0.5 * int.structure.atom_list[ii].mass  * dot(vel_half[ii,:],vel_half[ii,:])
		end
		T = 2*KE / (dof * int.units.Kb)
		zeta += 0.5 * int.delta_t * (T/T_tar - 1.0)/tau^2
    	#Update velocity
    	for ii in 1:int.structure.noa
    		int.structure.atom_list[ii].vel = (vel_half[ii,:] + 0.5 * int.structure.atom_list[ii].force/int.structure.atom_list[ii].mass* int.units.force_factor * int.delta_t )/(1 + 0.5 * int.delta_t * zeta )
		end
		int.temp_vec[i] = T
		int.zeta_vec[i] = zeta
		put_atoms_back_to_box(int.structure)
		if int.output_filename != "None"
			if i%int.steps_save == 0
				Jumol.write_box(int.structure.reader,int.structure,i)
			end
		end
	end
	## close output file
	if int.output_filename != "None"
		close_file(int.structure.reader)
	end
	print("\n")
end


function run_npt(int::Integ,num_steps)
	println("... npt not yet implemented ...")
end


function get_temp_vec(int::Integ,t_beg,t_end,num_steps)
	temp_target_vec = zeros(num_steps)
	for i in 1:num_steps
		temp_target_vec[i] = t_beg + (t_end-t_beg)/float(num_steps) * (i-1)
	end
	return temp_target_vec
end

## Bolzmann distribution
function set_atom_velocities(int::Integ, order_of_problem, T)
	for atom in int.structure.atom_list
		mu_vel = 0.0
		sigma_vel = sqrt(int.units.Kb * T/atom.mass)
		dist = Normal(mu_vel,sigma_vel)
		if order_of_problem == 2
			atom.vel = [rand(dist), rand(dist), 0.0]
		elseif order_of_problem == 3
			atom.vel = [rand(dist), rand(dist), rand(dist)]
		end
	end
end
