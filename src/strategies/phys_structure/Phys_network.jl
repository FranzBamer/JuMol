## Generate and modify the physical network
######################################
######################################
######################################
# physical network structure
#
#
#
#
#
# (c) Franz Bamer March-2022
######################################

include("../../structure/Structure.jl")
include("Node.jl")
include("Bond.jl")
include("../../calc/min.jl")


## constructor of the physical network (empty at the beginning)
mutable struct Phys_network
    dn # dn ... dual_network
    Phys_network(dn) = new(dn)
    node_list
    bond_list
    ring_list
	##
	physstruc
	physcalc
	##
	silicastruc
	silicacalc
	##
	include_lammps::Bool
	lammps_executable::String
	sig::Int64
end

function dual2phys(pn::Phys_network)
    ## build network structure from the dual
    pn.node_list = Vector()
    for i in 1:length(pn.dn.dual_triplet_list)
        tr = pn.dn.dual_triplet_list[i]
        num1 = tr.node1.number
        num2 = tr.node2.number
        num3 = tr.node3.number
        #
        c1 = pn.dn.dualstruc.atom_list[num1].pos[1:2]
        c2 = c1 - get_distance_between2atoms(pn.dn.dualstruc,num1,num2)[1:2]
        c3 = c2 - get_distance_between2atoms(pn.dn.dualstruc,num2,num3)[1:2]
        #
        center_coord = (c1+c2+c3)*1.0/3.0
        push!(pn.node_list,Node(i,center_coord[1],center_coord[2],1))
        find_neigh_triplets(tr,pn.dn.dual_triplet_list)
    end
    pn.bond_list = Vector()
    bond_cntr = 1
    for i in 1:length(pn.dn.dual_triplet_list)
        ctr = pn.dn.dual_triplet_list[i]
        for ii in 1:length(ctr.neigh_triplet_list)
            num_neigh_tr = pn.dn.dual_triplet_list[i].neigh_triplet_list[ii]
            if num_neigh_tr > ctr.number
                bond = Bond(bond_cntr,pn.node_list[ctr.number],pn.node_list[num_neigh_tr])
                push!(pn.bond_list,bond)
                bond_cntr += 1
            end
        end
    end
	## disturb one atom
	#pn.node_list[1].pos_x += 1.0
	#pn.node_list[1].pos_y += 1.0
end

## build molstruc of the network (harmonic potential)
function build_phys_molstruc(pn::Phys_network)
    ##
	pn.physstruc = Structure()
    pn.physstruc.rc = pn.dn.delta_x*6.0 # so that only the nearest neighbor is in the list
    pn.physstruc.pbx = 1
    pn.physstruc.pby = 1
	pn.physstruc.linked_cells_bool = false
    initialize_structure_objects(pn.physstruc)
	# load nodes into the structure
	for i in 1:length(pn.node_list)
		n = pn.node_list[i]
		add_atom_by_hand(pn.physstruc,i,n.type,n.pos_x,n.pos_y,0.0,0.0,0.0,0.0)
	end
	# create the box
	create_box_by_hand(pn.physstruc,pn.dn.lx,pn.dn.ly,1.0,0.0,0.0,0.0)
	## initialize the dual-molstructure
	initialize_structure(pn.physstruc)
	## manipulate bond structure to the network neighbors
	for i in 1:length(pn.dn.dual_triplet_list)
		atom = pn.physstruc.atom_list[i]
		atom.neighbor_indices = pn.dn.dual_triplet_list[i].neigh_triplet_list
	end
	## put atoms back to box
	put_atoms_back_to_box(pn.physstruc)
	update_node_pos(pn)
	## molcalc
	pn.physcalc = Calc(pn.physstruc)
	## harmonic ring-ring potential
	initialize_potential(pn.physcalc,7)
	pn.physcalc.potential.bl = pn.dn.bond_length
	set_ring_ring_params(pn.physcalc.potential)
	flush(stdout)
end

## minimize physical network (harmonic)
function minimize_phys_network(pn::Phys_network,num_steps=2000)
	Min_harm = Min(pn.physstruc,pn.physcalc)
	Min_harm.alpha_min = 1.0e-2
	#run_cg(Min_harm,1000,1.0e-5,[0],false)
	run_sd(Min_harm,1.0e-3,num_steps,[0],false) # full update -> false, because the the bonds have to be defined
end

## update positions from physstruc
function update_node_pos(pn::Phys_network)
	for i in 1:length(pn.node_list)
		pn.node_list[i].pos_x = pn.physstruc.atom_list[i].pos[1]
		pn.node_list[i].pos_y = pn.physstruc.atom_list[i].pos[2]
	end
end

## calculate a ring center
function calculate_ring_center(pn::Phys_network, adjacent_triplet_list)
	center_pos = zeros(2)
	num1 = adjacent_triplet_list[1]
	pos = pn.physstruc.atom_list[num1].pos[1:2]
	center_pos += pos
	for i in 2:length(adjacent_triplet_list)
		num1 = adjacent_triplet_list[i-1]
		num2 = adjacent_triplet_list[i]
		pos = pos - get_distance_between2atoms(pn.physstruc,num1,num2)[1:2]
		center_pos += pos
	end
	return center_pos * 1.0/float(length(adjacent_triplet_list))
end

## build molstruc of the silica network (silica potential)
function build_silica_molstruc(pn::Phys_network)
    ##
	pn.silicastruc = Structure()
    pn.silicastruc.rc = 10.17
    pn.silicastruc.pbx = 1
    pn.silicastruc.pby = 1
    initialize_structure_objects(pn.silicastruc)
	# load nodes into the structure
	atom_cntr = 1
	for i in 1:length(pn.node_list)
		n = pn.node_list[i]
		add_atom_by_hand(pn.silicastruc,atom_cntr,1,n.pos_x,n.pos_y,0.0,0.0,0.0,0.0)
		atom_cntr += 1
	end
	for i in 1:length(pn.bond_list)
		bond = pn.bond_list[i]
		n1 = bond.node1.number
		n2 = bond.node2.number
		pos1 = pn.physstruc.atom_list[n1].pos[1:2]
		pos2 = pos1 - get_distance_between2atoms(pn.physstruc,n1,n2)[1:2]
		center_coord = (pos1 + pos2) * 0.5
		add_atom_by_hand(pn.silicastruc,atom_cntr,2,
		                 center_coord[1],center_coord[2],0.0,0.0,0.0,0.0)
		atom_cntr += 1
	end
	# create the box
	create_box_by_hand(pn.silicastruc,pn.dn.lx,pn.dn.ly,1.0,0.0,0.0,0.0)
	## initialize the dual-molstructure
	initialize_structure(pn.silicastruc)
	## put atoms back to box
	put_atoms_back_to_box(pn.silicastruc)
	## molcalc
	pn.silicacalc = Calc(pn.silicastruc)
	## harmonic ring-ring potential
	initialize_potential(pn.silicacalc,13)
	flush(stdout)
end

## minimize physical network (harmonic)
function minimize_silica_network(pn::Phys_network, acc_factor=1.001, num_steps=2000, tolerance=1.0e-5)
	println(" ")
	Min_silica = Min(pn.silicastruc,pn.silicacalc)
	Min_silica.alpha_min = 1.0e-1
    Min_silica.acc_factor = acc_factor
	run_cg(Min_silica, num_steps, tolerance, [0], true)
end


## update positions from silicastruc
function update_silica_node_pos(pn::Phys_network)
	for i in 1:length(pn.node_list)
		pn.node_list[i].pos_x = pn.silicastruc.atom_list[i].pos[1]
		pn.node_list[i].pos_y = pn.silicastruc.atom_list[i].pos[2]
	end
end

## consistency check for the silica network: Si-Si -> type_neighbor 1, Si-O -> type neighbor 2
function silica_consistency_check(pn::Phys_network,type_neighbor=2,dist_neigh=1.8)
	check = true
	for i in 1:length(pn.node_list)
		center_atom = pn.silicastruc.atom_list[i]
		if center_atom.type == 1
			num_neighbors = 0
			for ii in 1:length(center_atom.neighbor_indices)
				num_neigh = center_atom.neighbor_indices[ii]
				type_neigh = pn.silicastruc.atom_list[num_neigh].type
				if type_neigh==type_neighbor
					dist = get_distance_between2atoms(pn.silicastruc,center_atom.number,num_neigh)[4]
					if dist < dist_neigh
						num_neighbors += 1
					end
				end
			end
			if num_neighbors != 3
				check =false
				break
			end
		end
	end
	return check
end


#### function that runs lammps for the minimization
## minimize silica structure in LAMMPS
function minimize_silica_network_lammps_plus_coord_check(pn::Phys_network)
	## write temporary input-file
	open_file(pn.silicastruc.reader, "sample_tmp_"*string(pn.sig)*".twodsilica", "w")
	write_lmp_box(pn.silicastruc.reader, pn.silicastruc, 1)
	close_file(pn.silicastruc.reader)
	## check the minimum atomic distance here:
	# if it is smaller than 0.6 (smallest table number) then coordination_check
	# is false
	min_dist = find_minimum_pair_distance(pn,pn.silicastruc)
	if min_dist < 0.61
		coordination_check = false
	else
		## run lammps minimization
		lmp_serial = pn.lammps_executable
		lmp_input_file = "relax_"*string(pn.sig)*".in"
		cmd_lammps_min = `$lmp_serial -in $lmp_input_file -log none`
		run(cmd_lammps_min)
		# go back into the original calculation folder
		#cd("../")
		## load lammps result back into JuMol
		pn.silicastruc = Structure()
		pn.silicastruc.rc = 10.0
		pn.silicastruc.rskin = 0.5
		pn.silicastruc.pbx = 1
		pn.silicastruc.pby = 1
		pn.silicastruc.linked_cells_bool = true
		initialize_structure_objects(pn.silicastruc)
		lammps_output_file = "rel_sample_"*string(pn.sig)*".lammpstrj"
		read_lammpstrj(pn.silicastruc,lammps_output_file)
		initialize_structure(pn.silicastruc)
		pn.silicacalc = Calc(pn.silicastruc)
		initialize_potential(pn.silicacalc,13)
		## check the coordination number
		coordination_check = silica_consistency_check(pn)
		#readline()
		## clean-up
		input_file_silica = "sample_tmp_"*string(pn.sig)*".twodsilica"
		run(`rm $input_file_silica`)
		output_file_silica = "rel_sample_"*string(pn.sig)*".lammpstrj"
		run(`rm $output_file_silica`)
		coordination_file = "coordination_"*string(pn.sig)*".dat"
		run(`rm $coordination_file`)
		log_lmp_file = "log_"*string(pn.sig)*".lammps"
		run(`rm $log_lmp_file`)
		#readline()
	end
	return coordination_check
end


#### function that uses a lammps output file as information
## consitency check
function silica_consistency_check_lammps(pn::Phys_network)
	check = true
	file = open("coordination_"*string(pn.sig)*".dat","r")
	lines = readlines(file)
	close(file)
	## do not read head (10 lines)
	cntr_line = 10
	while cntr_line < length(lines)
		line_vec = readdlm(IOBuffer(reader.lines[cntr_line]))
		num_atom = Int(line_vec[1])
		type_atom = Int(line_vec[2])
		coordination = Int(line_vec[3])
		# Si atom must have 3 neighbors
		if type_atom == 1
			if coordination != 3
				check = false
			end
		end
		# O atom must have 2 neighbors
		if type_atom == 2
			if coordination != 2
				check = false
			end
		end
		#
		cntr_line += 1
	end
	return check
end

## finding the minimum distance between any possible atom pair in the sample
function find_minimum_pair_distance(pn::Phys_network,silicastruc)
	min_dist = 999.9
	for i in 1:length(silicastruc.atom_list)
		center_atom = silicastruc.atom_list[i]
		for ii in center_atom.neighbor_indices
			dist = get_distance_between2atoms(silicastruc,i,ii)[4]
			if dist < min_dist
				min_dist = dist
			end
			# println(" ")
			# println("cenat: ", i)
			# println("neighat: ", ii)
			# println(dist)
		end
	end
	# println("###################################")
	# println("min_dist: ", min_dist)
	#readline()
	return min_dist
end
