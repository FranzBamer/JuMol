## Distance calculations
######################################
######################################
######################################
# full distance updates
# verlet distance updates
# linked cell distance updates
#
# (c) Franz Bamer, Ivy Wu Mai-2022
######################################


## calculate all distances for the whole model
function update_distances(structure::Structure)
	for i in 1:structure.noa
		if structure.linked_cells_bool == true
			## calculate distance update using the linked cell lists
			for cell_id in structure.cell_neighbour_array[structure.atom_list[i].atom_cell_index_vec[4],:]
				if cell_id != 0 && length(structure.cell_list[cell_id].atoms_in_cell_list)>0
					for ii in structure.cell_list[cell_id].atoms_in_cell_list
						if i != ii
							dist_vec = calc_distance_betw2atoms(structure,min(i,ii),max(i,ii))
							structure.atom_list[min(i,ii)].distances[:,max(i,ii)-min(i,ii)] = dist_vec
						end
					end
				end
			end
		else
			## calculate full distance update
			for ii in i+1:structure.noa
				dist_vec = calc_distance_betw2atoms(structure, i, ii)
				structure.atom_list[i].distances[:,ii-i] = dist_vec
			end
		end
	end
end

## update neighbor lists for every atom
function update_neighbor_lists(structure::Structure)
	for i in 1:structure.noa
		## erase old list
		structure.atom_list[i].neighbor_indices = Vector{Int64}()
		if structure.linked_cells_bool == true
			## neighbor list update considering linked cell list
			for cell_id in structure.cell_neighbour_array[structure.atom_list[i].atom_cell_index_vec[4],:]
				if cell_id != 0 && length(structure.cell_list[cell_id].atoms_in_cell_list)>0
					for ii in structure.cell_list[cell_id].atoms_in_cell_list
						if i != ii
							if get_distance_between2atoms(structure,min(i,ii),max(i,ii))[4] < structure.rc+structure.rskin
								push!(structure.atom_list[i].neighbor_indices,ii)
							end
						end
					end
				end
			end
		else
			## neigbor list update full
			for ii in 1:structure.noa
				if i != ii
					if get_distance_between2atoms(structure,i,ii)[4] < structure.rc+structure.rskin
						push!(structure.atom_list[i].neighbor_indices,ii)
					end
				end
			end
		end
	end
end

## update distances in within the cutoff radius
function update_distances_partly(structure::Structure)
	for i in 1:structure.noa
		update_distances_partly_atom(structure,i)
	end
end

## update the distances for one atom within its neighor list
function update_distances_partly_atom(structure::Structure,atom_num)
	for i in structure.atom_list[atom_num].neighbor_indices
		if i > atom_num
			dist_vec = calc_distance_betw2atoms(structure, atom_num, i)
			structure.atom_list[atom_num].distances[:,i-atom_num] = dist_vec
		end
	end
end

## update the cell list and information
function update_cell(structure::Structure)
	# modify the cutoff radius for shear condition
	tan_gamma = structure.box.lxy/structure.box.ly
	cos_gamma = sqrt(structure.box.lxy * structure.box.lxy + structure.box.ly * structure.box.ly) / structure.box.ly
	r_cutoff = (structure.rc+structure.rskin)/cos_gamma
	# subdivide global box into n smaller cells
	n_cell_x = floor(structure.box.lx/r_cutoff)
	n_cell_y = floor(structure.box.ly/r_cutoff)
	n_cell_z = floor(structure.box.lz/r_cutoff)
	# 2D case
	if n_cell_x == 0
		n_cell_x = 1
	end
	if n_cell_y == 0
		n_cell_y = 1
	end
	if n_cell_z == 0
		n_cell_z = 1
	end
	# cell length
	l_cell_x = structure.box.lx/n_cell_x
	l_cell_y = structure.box.ly/n_cell_y
	l_cell_z = structure.box.lz/n_cell_z
	# total number of cells
	n_cell_xyz = n_cell_x*n_cell_y*n_cell_z
	# erase the original cell list
	structure.cell_list = Vector()
	#
	for i in 0:n_cell_x-1
		for ii in 0:n_cell_y-1
			for iii in 0:n_cell_z-1
				cell_id = Int(i*n_cell_y*n_cell_z + ii*n_cell_z + iii + 1)
				push!(structure.cell_list, Cell(cell_id))
				# erase the original cell atom list
				structure.cell_list[cell_id].atoms_in_cell_list = Vector{Int64}()
				#println(cell_id)
			end
		end
	end
end

## create list of linked cells
function linked_cell_list(structure::Structure)
	# modify the cutoff radius for shear condition
	tan_gamma = structure.box.lxy/structure.box.ly
	cos_gamma = sqrt(structure.box.lxy * structure.box.lxy + structure.box.ly * structure.box.ly) / structure.box.ly
	r_cutoff = (structure.rc+structure.rskin)/cos_gamma
	# subdivide global box into n smaller cells
	n_cell_x = floor(structure.box.lx/r_cutoff)
	n_cell_y = floor(structure.box.ly/r_cutoff)
	n_cell_z = floor(structure.box.lz/r_cutoff)
	# 2D case
	if n_cell_x == 0
		n_cell_x = 1
	end
	if n_cell_y == 0
		n_cell_y = 1
	end
	if n_cell_z == 0
		n_cell_z = 1
	end
	# cell length  !Fixme! if 2d, there will be a zero in n_cell_z
	l_cell_x = structure.box.lx/n_cell_x
	l_cell_y = structure.box.ly/n_cell_y
	l_cell_z = structure.box.lz/n_cell_z
	# total number of cells
	n_cell_xyz = n_cell_x*n_cell_y*n_cell_z
	#
	# assign every atom to cell and link cell index with atom index
	for i in 1:structure.noa
		# put atoms into original box for shear
		atom_list_x = structure.atom_list[i].pos[1] - tan_gamma * structure.atom_list[i].pos[2]
		atom_list_y = structure.atom_list[i].pos[2]
		atom_list_z = structure.atom_list[i].pos[3]
		#
		# ckeck to which cell index atom i belongs
		structure.atom_list[i].atom_cell_index_vec = zeros(4)
		structure.atom_list[i].atom_cell_index_vec[1] = Int(floor(atom_list_x/l_cell_x))
		structure.atom_list[i].atom_cell_index_vec[2] = Int(floor(atom_list_y/l_cell_y))
		structure.atom_list[i].atom_cell_index_vec[3] = Int(floor(atom_list_z/l_cell_z))
		#
		# boundary correction
		if structure.atom_list[i].atom_cell_index_vec[1]<0.0
			structure.atom_list[i].atom_cell_index_vec[1] = n_cell_x-1
		elseif structure.atom_list[i].atom_cell_index_vec[1]>=n_cell_x
			structure.atom_list[i].atom_cell_index_vec[1] = 0
		end
		if structure.atom_list[i].atom_cell_index_vec[2]<0.0
			structure.atom_list[i].atom_cell_index_vec[2] = n_cell_y-1
		elseif structure.atom_list[i].atom_cell_index_vec[2]>=n_cell_y
			structure.atom_list[i].atom_cell_index_vec[2] = 0
		end
		if structure.atom_list[i].atom_cell_index_vec[3]<0.0
			structure.atom_list[i].atom_cell_index_vec[3] = n_cell_z-1
		elseif structure.atom_list[i].atom_cell_index_vec[3]>=n_cell_z
			structure.atom_list[i].atom_cell_index_vec[3] = 0
		end
		#
		# convert cell index vec into scalar
		cell_index = structure.atom_list[i].atom_cell_index_vec[1]*n_cell_y*n_cell_z + structure.atom_list[i].atom_cell_index_vec[2]*n_cell_z + structure.atom_list[i].atom_cell_index_vec[3] + 1
		structure.atom_list[i].atom_cell_index_vec[4] = Int(cell_index)
		#
		# the last one goes to the header
		push!(structure.cell_list[Int(cell_index)].atoms_in_cell_list, Int(i))
	end
end

## create list of neighboring cells
function construct_neigh_cell_list(structure::Structure)
	# modify the cutoff radius for shear condition
	tan_gamma = structure.box.lxy/structure.box.ly
	cos_gamma = sqrt(structure.box.lxy * structure.box.lxy + structure.box.ly * structure.box.ly) / structure.box.ly
	r_cutoff = (structure.rc+structure.rskin)/cos_gamma
	# subdivide global box into n smaller cells
	n_cell_x = floor(structure.box.lx/r_cutoff)
	n_cell_y = floor(structure.box.ly/r_cutoff)
	n_cell_z = floor(structure.box.lz/r_cutoff)
	# 2D case
	if n_cell_x == 0
		n_cell_x = 1
	end
	if n_cell_y == 0
		n_cell_y = 1
	end
	if n_cell_z == 0
		n_cell_z = 1
	end
	# cell length
	l_cell_x = structure.box.lx/n_cell_x
	l_cell_y = structure.box.ly/n_cell_y
	l_cell_z = structure.box.lz/n_cell_z
	#total number of cells
	n_cell_xyz = n_cell_x*n_cell_y*n_cell_z
	#
    # construct lists of neighbouring cells
	neighboring_cell_array = zeros(Int64,Int(n_cell_xyz),27)
	structure.cell_neighbour_array = zeros(Int64,Int(n_cell_xyz),27)
	# current cell
	for i in 0:n_cell_x-1
		for ii in 0:n_cell_y-1
			for iii in 0:n_cell_z-1
				c = Int(i*n_cell_y*n_cell_z + ii*n_cell_z + iii + 1)
				# neighbour cells
				count_neigh = 1
				for j in i-1:i+1
					for jj in ii-1:ii+1
						for jjj in iii-1:iii+1
							#
							if structure.pbx == 1
								if j<0
									j = j + n_cell_x
								elseif j==n_cell_x
									j = 0
								end
							end
							#
							if structure.pby == 1
								if jj<0
									jj = jj + n_cell_y
								elseif jj==n_cell_y
									jj = 0
								end
							end
							#
							if structure.pbz == 1
								if jjj<0
									jjj = jjj + n_cell_z
								elseif jjj==n_cell_z
									jjj = 0
								end
							end
							#
							if j>=0 && jj>=0 && jjj>=0 && j<n_cell_x && jj<n_cell_y && jjj<n_cell_z
								c1 = Int(j*n_cell_y*n_cell_z + jj*n_cell_z + jjj +1)
								neighboring_cell_array[c,count_neigh] = c1
								structure.cell_neighbour_array[c,count_neigh] = c1
							end
							count_neigh += 1
							#
						end
					end
				end
			end
		end
	end
end

## calculate the distance between two atoms
function calc_distance_betw2atoms(structure::Structure,num1,num2)
	pos1 = structure.atom_list[num1].pos
	pos2 = copy(structure.atom_list[num2].pos)
	dist_vec = zeros(4)
	dist_vec[1:3] = pos2 - pos1
	## periodic boundary conditions in e1
	if structure.pbx == 1
		if dot(dist_vec[1:3],structure.box.e1) > structure.box.l1*0.5
			pos2 = pos2 - structure.box.h1
		end
		if dot(dist_vec[1:3],structure.box.e1) < -structure.box.l1*0.5
			pos2 = pos2 + structure.box.h1
		end
	end
	## periodic boundary conditions in e2
	if structure.pby == 1
		if dot(dist_vec[1:3],structure.box.e2) > structure.box.l2*0.5
			pos2 = pos2 - structure.box.h2
		end
		if dot(dist_vec[1:3],structure.box.e2) < -structure.box.l2*0.5
			pos2 = pos2 + structure.box.h2
		end
	end
	## periodic boundary conditions in e3
	if structure.pbz == 1
		if dot(dist_vec[1:3],structure.box.e3) > structure.box.l3*0.5
			pos2 = pos2 - structure.box.h3
		end
		if dot(dist_vec[1:3],structure.box.e3) < -structure.box.l3*0.5
			pos2 = pos2 + structure.box.h3
		end
	end
	if num1 < num2
		dist_vec[1:3] = pos1 - pos2
	else
		dist_vec[1:3] = pos2 - pos1
	end
	dist_vec[4] = norm(dist_vec[1:3])
	return dist_vec
end

## get the distance between two atoms
function get_distance_between2atoms(structure::Structure, num1, num2)
	dist_vec = zeros(4)
	if num1 <= num2
		dist_vec = structure.atom_list[num1].distances[:,num2-num1]
	else
		dist_vec = structure.atom_list[num2].distances[:,num1-num2]
		dist_vec[1:3] = dist_vec[1:3]*(-1.0)
	end
	return dist_vec
end
