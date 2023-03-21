## Generate and modify a dual network
######################################
######################################
######################################
# dual network structure
#
#
#
#
#
# (c) Franz Bamer Nov-2020
######################################

include("Dual_node.jl")
include("Dual_bond.jl")
include("Dual_triplet.jl")
include("../../structure/Structure.jl")
include("../../calc/calc.jl")
include("../../calc/min.jl")

## constructor of the dual lattice (num_y must be even)
mutable struct Dual_network
    bond_length::Float64
    Dual_network(bond_length) = new(bond_length)
	num_x::Int64
    num_y::Int64
	## lattice constants
	a::Float64
	delta_x::Float64
	delta_y::Float64
	## dual box
	lx::Float64
	ly::Float64
	lz::Float64
	## dual node list
	dual_node_list
	dual_bond_list
	dual_triplet_list
	## molstruc and molcalc of the dual network
	dualstruc
	dualcalc
	dualmin
	##
	r_list
end



## start with the hexagonal lattice
function build_starting_dual_hex(dn::Dual_network,num_x,num_y,r_list=[4,5,6,7,8,9,10])
	## general allowed ring_size_list
	dn.r_list = r_list
	##
	dn.num_x = num_x
	dn.num_y = num_y
	# lattice parameters
	dn.a = 2.0*dn.bond_length*cos(pi/180.0*30.0)
	dn.delta_x = dn.a
	dn.delta_y = dn.a*sin(60.0*pi/180.0)
	# size of the dual box
	dn.lx = dn.num_x*dn.delta_x
	dn.ly = dn.num_y*dn.delta_y
	dn.lz = 1.0
	println(" - building initial hexagonal network - ")
	## creating dual node list of the structure
	dn.dual_node_list = Vector()
	cntr = 1
	for i in 1:dn.num_y
		pos_y = (i-1)*dn.delta_y
		pos_x_start = 0.0 + float((i-1)%2)*dn.a*0.5
		for ii in 1:dn.num_x
			## add triangle dual node
			push!(dn.dual_node_list,Dual_node(cntr, 6, pos_x_start+(ii-1)*dn.delta_x, pos_y))
			cntr += 1
		end
	end
	## create dual bonds and dual triplets of the hexagonal structure
	dn.dual_bond_list = Vector()
	dn.dual_triplet_list = Vector()
	## first set of triplets
	bond_cntr = 1
	triplet_cntr = 1
	## lying triangles
	for j in 1:dn.num_y-1
		for i in 1:dn.num_x
			if j%2 != 0
				# define node numbers
				node1 = dn.dual_node_list[dn.num_x*(j-1)+i]
				if i<dn.num_x
					node2 = dn.dual_node_list[dn.num_x*(j-1)+i+1]
				else
					node2 = dn.dual_node_list[dn.num_x*(j-1)+1]
				end
				node3 = dn.dual_node_list[dn.num_x*j+ i]
				## bond1 1 in triplet
				bond1 = Dual_bond(bond_cntr,node1,node2)
				bond_cntr += 1
				## bond 2 in triplet
				bond2 = Dual_bond(bond_cntr,node2,node3)
				bond_cntr += 1
				## bond 3 in triplet
				bond3 = Dual_bond(bond_cntr,node3,node1)
				bond_cntr += 1
				##
				push!(dn.dual_bond_list,bond1)
				push!(dn.dual_bond_list,bond2)
				push!(dn.dual_bond_list,bond3)
				## triplet
				triplet = Dual_triplet(triplet_cntr,bond1,bond2,bond3)
				find_triangle_nodes(triplet)
				push!(dn.dual_triplet_list,triplet)
				triplet_cntr += 1
			else # j is an even number
				## define node numbers
				if i==1
					node1 = dn.dual_node_list[dn.num_x*(j-1)+dn.num_x]
				else
					node1 = dn.dual_node_list[dn.num_x*(j-1)+i-1]
				end
				node2 = dn.dual_node_list[dn.num_x*(j-1)+i]
				node3 = dn.dual_node_list[dn.num_x*j+ i]
				## bond1 1 in triplet
				bond1 = Dual_bond(bond_cntr,node1,node2)
				bond_cntr += 1
				## bond 2 in triplet
				bond2 = Dual_bond(bond_cntr,node2,node3)
				bond_cntr += 1
				## bond 3 in triplet
				bond3 = Dual_bond(bond_cntr,node3,node1)
				bond_cntr += 1
				##
				push!(dn.dual_bond_list,bond1)
				push!(dn.dual_bond_list,bond2)
				push!(dn.dual_bond_list,bond3)
				## triplet
				triplet = Dual_triplet(triplet_cntr,bond1,bond2,bond3)
				find_triangle_nodes(triplet)
				push!(dn.dual_triplet_list,triplet)
				triplet_cntr += 1
			end
		end
	end
	## standing triangles
	for j in 1:dn.num_y-1
		for i in 1:dn.num_x
			# even rows
			if j%2 != 0 && j<dn.num_y-1
				##
				bond1 = dn.dual_bond_list[dn.num_x*3*(j-1)+3*(i)]
				bond2 = dn.dual_bond_list[dn.num_x*3*j+3*(i-1)+1]
				if i==1
					bond3 = dn.dual_bond_list[dn.num_x*3*j-1]
				else
					bond3 = dn.dual_bond_list[dn.num_x*3*(j-1)+3*(i-1)-1]
				end
				## triplet
				triplet = Dual_triplet(triplet_cntr,bond1,bond2,bond3)
				find_triangle_nodes(triplet)
				push!(dn.dual_triplet_list,triplet)
				triplet_cntr += 1
			end
			# odd rows
			if j%2 == 0 && j<dn.num_y-1
				##
				if i == dn.num_x
					bond1 = dn.dual_bond_list[dn.num_x*3*(j-1)+3]
					bond2 = dn.dual_bond_list[dn.num_x*3*j+i*3-2]
					bond3 = dn.dual_bond_list[dn.num_x*3*(j-1)+i*3-1]
				else
					bond1 = dn.dual_bond_list[dn.num_x*3*(j-1)+(i+1)*3]
					bond2 = dn.dual_bond_list[dn.num_x*3*j+i*3-2]
					bond3 = dn.dual_bond_list[dn.num_x*3*(j-1)+i*3-1]
				end
				## triplet
				triplet = Dual_triplet(triplet_cntr,bond1,bond2,bond3)
				find_triangle_nodes(triplet)
				push!(dn.dual_triplet_list,triplet)
				triplet_cntr += 1
			end
		end
	end
	## row of standing triangles for row number self.num_y-1
	for i in 1:dn.num_x
		if i==1
			bond1 = dn.dual_bond_list[dn.num_x*3*(dn.num_y-2)+3*i]
			node1 = dn.dual_node_list[dn.num_x*(dn.num_y-1)+i]
			node2 = dn.dual_node_list[dn.num_x*dn.num_y]
			bond2 = Dual_bond(dn.num_x*3*(dn.num_y-1)+1,node1,node2)
			push!(dn.dual_bond_list,bond2)
			bond3 = dn.dual_bond_list[dn.num_x*3*(dn.num_y-1)-1]
		else
			bond1 = dn.dual_bond_list[dn.num_x*3*(dn.num_y-2)+3*i]
			node1 = dn.dual_node_list[dn.num_x*(dn.num_y-1) + i]
			node2 = dn.dual_node_list[dn.num_x*(dn.num_y-1) + i - 1]
			bond2 = Dual_bond(dn.num_x*3*(dn.num_y-1)+i,node1,node2)
			push!(dn.dual_bond_list,bond2)
			bond3 = dn.dual_bond_list[dn.num_x*3*(dn.num_y-2)+3*(i-1)-1]
		end
		bond_cntr += 1
		## triplet
		triplet = Dual_triplet(triplet_cntr,bond1,bond2,bond3)
		find_triangle_nodes(triplet)
		push!(dn.dual_triplet_list,triplet)
		triplet_cntr += 1
	end
	## last row lying triangles
	for i in 1:dn.num_x
		if i<dn.num_x
			node1 = dn.dual_node_list[dn.num_x*(dn.num_y-1)+i]
			node2 = dn.dual_node_list[dn.num_x*(dn.num_y-1)+i+1]
			node3 = dn.dual_node_list[i+1]
			bond1 = dn.dual_bond_list[3*dn.num_x*(dn.num_y-1)+i+1]
			bond2 = Dual_bond(bond_cntr,node2,node3)
			push!(dn.dual_bond_list,bond2)
			bond_cntr += 1
			bond3 = Dual_bond(bond_cntr,node3,node1)
			push!(dn.dual_bond_list,bond3)
			bond_cntr += 1
		else
			node1 = dn.dual_node_list[dn.num_x*(dn.num_y-1)+i]
			node2 = dn.dual_node_list[dn.num_x*(dn.num_y-1)+1]
			node3 = dn.dual_node_list[1]
			bond1 = dn.dual_bond_list[3*dn.num_x*(dn.num_y-1)+1]
			bond2 = Dual_bond(bond_cntr,node2,node3)
			push!(dn.dual_bond_list,bond2)
			bond_cntr += 1
			bond3 = Dual_bond(bond_cntr,node3,node1)
			push!(dn.dual_bond_list,bond3)
			bond_cntr += 1
		end
		## triplet
		triplet = Dual_triplet(triplet_cntr,bond1,bond2,bond3)
		find_triangle_nodes(triplet)
		push!(dn.dual_triplet_list,triplet)
		triplet_cntr += 1
	end
	## last row standing triangles
	for i in 1:dn.num_x
		if i==1
			bond1 = dn.dual_bond_list[3*dn.num_x*(dn.num_y-1)+dn.num_x+2]
			bond2 = dn.dual_bond_list[1]
			bond3 = dn.dual_bond_list[3*dn.num_x*dn.num_y-1]
		else
			bond1 = dn.dual_bond_list[3*dn.num_x*(dn.num_y-1)+dn.num_x+2*i]
			bond2 = dn.dual_bond_list[i*3-2]
			bond3 = dn.dual_bond_list[3*dn.num_x*(dn.num_y-1)+dn.num_x+2*i-3]
		end
		## triplet
		triplet = Dual_triplet(triplet_cntr,bond1,bond2,bond3)
		find_triangle_nodes(triplet)
		push!(dn.dual_triplet_list,triplet)
		triplet_cntr += 1
	end
	##
	for i in 1:length(dn.dual_node_list)
		get_adjacent_neigh_nodes(dn.dual_node_list[i],dn.dual_node_list,dn.dual_bond_list)
	end
	##
	build_dual_molstruc(dn)
	## find the neighboring triangles
	for i in 1:length(dn.dual_triplet_list)
		find_neigh_triplets(dn.dual_triplet_list[i],dn.dual_triplet_list)
	end
	## ring calculation
	## calculate neighboring triplets to every node
	for i in 1:length(dn.dual_node_list)
		find_adjacent_neigh_triplets_ordered(dn.dual_node_list[i],dn.dual_triplet_list)
	end
end



## build the dual molstruc
function build_dual_molstruc(dn::Dual_network)
	##
	dn.dualstruc = Structure()
    dn.dualstruc.rc = dn.delta_x*3.0 # so that only the nearest neighbor is in the list
    dn.dualstruc.pbx = 1
    dn.dualstruc.pby = 1
	dn.dualstruc.linked_cells_bool = false
    initialize_structure_objects(dn.dualstruc)
	# load nodes into the structure
	for i in 1:length(dn.dual_node_list)
		n = dn.dual_node_list[i]
		add_atom_by_hand(dn.dualstruc,i,n.type,n.pos_x,n.pos_y,0.0,0.0,0.0,0.0)
	end
	# create the box
	create_box_by_hand(dn.dualstruc,dn.lx,dn.ly,1.0,0.0,0.0,0.0)
	## initialize the dual-molstructure
	initialize_structure(dn.dualstruc)
	## manipulate bond structure to the network neighbors
	for i in 1:length(dn.dual_node_list)
		atom = dn.dualstruc.atom_list[i]
		node = dn.dual_node_list[i]
		atom.neighbor_indices = node.adjacent_node_list
	end

	## molcalc
	dn.dualcalc = Calc(dn.dualstruc)
	## harmonic ring-ring potential
	initialize_potential(dn.dualcalc,7)
	dn.dualcalc.potential.bl = dn.bond_length
	set_ring_ring_params(dn.dualcalc.potential)
	flush(stdout)
end


## minimize the dual network
function minimize_dual_molstruc(dn::Dual_network,num_steps=2000)
	dn.dualmin = Min(dn.dualstruc,dn.dualcalc)
	dn.dualmin.alpha_min = 0.1
	run_sd(dn.dualmin,1.0e-3,num_steps,[0],false)
end

## update node coordinate information
function update_dual_network(dn::Dual_network)
	for i in 1:length(dn.dual_node_list)
		nd = dn.dual_node_list[i]
		nd.pos_x = dn.dualstruc.atom_list[i].pos[1]
		nd.pos_y = dn.dualstruc.atom_list[i].pos[2]
	end
end

## update dual node position
function phys2dual(dn::Dual_network,phys_network)
	for i in 1:length(dn.dual_node_list)
		dual_node = dn.dual_node_list[i]
		pos_center = calculate_ring_center(phys_network,dual_node.adjacent_triplet_list_ordered)
		dual_node.pos_x = pos_center[1]
		dual_node.pos_y = pos_center[2]
		dn.dualstruc.atom_list[i].pos[1:2] = pos_center
	end
	update_distances_partly(dn.dualstruc)
end

## switch a dual bond
function switch_bond(dn::Dual_network,bond_num)
	## define bond to switch
	bond = dn.dual_bond_list[bond_num]
	n1_old_num = bond.node1.number
	n2_old_num = bond.node2.number
	#println("nodes from which to switch:")
	#println(n1_old_num, ", ", n2_old_num)
	## find neighboring triangles
	adjacent_triplets = find_neighboring_triplets(bond,dn.dual_triplet_list)
	triplet1 = adjacent_triplets[1]
	triplet2 = adjacent_triplets[2]
	#println("corresponding triplet numbers:")
	#println(triplet1.number, ", ", triplet2.number)
	## find third node of the triangle 1, new number 1
	n1_new_num = find_third_node_of_triplet(triplet1,bond).number
	## find third node of the triangle 2, new number 2
	n2_new_num = find_third_node_of_triplet(triplet2,bond).number
	#println("switch to the new node numbers")
	#println(n1_new_num,", ",n2_new_num)
	## rotate the bond -> redefine new numbers
	bond.node1 = dn.dual_node_list[n1_new_num]
	bond.node2 = dn.dual_node_list[n2_new_num]
	## change the node types (ring sizes)
	dn.dual_node_list[n1_old_num].type -= 1
	dn.dual_node_list[n2_old_num].type -= 1
	dn.dual_node_list[n1_new_num].type += 1
	dn.dual_node_list[n2_new_num].type += 1
	## updating node neighbors
	get_adjacent_neigh_nodes(dn.dual_node_list[n1_old_num],dn.dual_node_list,dn.dual_bond_list)
	#println("dual node neighbors 1:")
	#println(dn.dual_node_list[n1_old_num].adjacent_node_list)
	get_adjacent_neigh_nodes(dn.dual_node_list[n2_old_num],dn.dual_node_list,dn.dual_bond_list)
	#println("dual node neighbors 2:")
	#println(dn.dual_node_list[n2_old_num].adjacent_node_list)
	get_adjacent_neigh_nodes(dn.dual_node_list[n1_new_num],dn.dual_node_list,dn.dual_bond_list)
	#println("dual node neighbors 3:")
	#println(dn.dual_node_list[n1_new_num].adjacent_node_list)
	get_adjacent_neigh_nodes(dn.dual_node_list[n2_new_num],dn.dual_node_list,dn.dual_bond_list)
	#println("dual node neighbors 4:")
	#println(dn.dual_node_list[n2_new_num].adjacent_node_list)
	## updating the triplets
	bond1 = find_bond_by_nodes(dn, dn.dual_node_list[n1_new_num],dn.dual_node_list[n2_new_num])
	# create new triplet 1
	bond21 = find_bond_by_nodes(dn, dn.dual_node_list[n1_new_num],dn.dual_node_list[n1_old_num])
	bond31 = find_bond_by_nodes(dn, dn.dual_node_list[n2_new_num],dn.dual_node_list[n1_old_num])
	dn.dual_triplet_list[triplet1.number] = Dual_triplet(triplet1.number,bond1,bond21,bond31)
	find_triangle_nodes(dn.dual_triplet_list[triplet1.number])
	# create new triplet 2
	bond22 = find_bond_by_nodes(dn, dn.dual_node_list[n1_new_num],dn.dual_node_list[n2_old_num])
	bond32 = find_bond_by_nodes(dn, dn.dual_node_list[n2_new_num],dn.dual_node_list[n2_old_num])
	dn.dual_triplet_list[triplet2.number] = Dual_triplet(triplet2.number,bond1,bond22,bond32)
	find_triangle_nodes(dn.dual_triplet_list[triplet2.number])
	# find neighbor triplets
	find_neigh_triplets(dn.dual_triplet_list[triplet1.number],dn.dual_triplet_list)
	#println("Neigh triplet list 1:")
	#println(dn.dual_triplet_list[triplet1.number].neigh_triplet_list)
	find_neigh_triplets(dn.dual_triplet_list[triplet2.number],dn.dual_triplet_list)
	#println("Neigh triplet list 2:")
	#println(dn.dual_triplet_list[triplet2.number].neigh_triplet_list)
	# update the neighboring triplets around the dual nodes that have changed
	find_adjacent_neigh_triplets_ordered(dn.dual_node_list[n1_old_num],dn.dual_triplet_list)
	#println("dual_ring_list 1:")
	#println(dn.dual_node_list[n1_old_num].adjacent_triplet_list_ordered)
	find_adjacent_neigh_triplets_ordered(dn.dual_node_list[n2_old_num],dn.dual_triplet_list)
	#println("dual_ring_list 2:")
	#println(dn.dual_node_list[n2_old_num].adjacent_triplet_list_ordered)
	find_adjacent_neigh_triplets_ordered(dn.dual_node_list[n1_new_num],dn.dual_triplet_list)
	#println("dual_ring_list 3:")
	#println(dn.dual_node_list[n1_new_num].adjacent_triplet_list_ordered)
	find_adjacent_neigh_triplets_ordered(dn.dual_node_list[n2_new_num],dn.dual_triplet_list)
	#println("dual_ring_list 4:")
	#println(dn.dual_node_list[n2_new_num].adjacent_triplet_list_ordered)
	## change dualstruc object
	# node1 (n1_old_num)
	atom = dn.dualstruc.atom_list[n1_old_num]
	node = dn.dual_node_list[n1_old_num]
	atom.type = node.type
	atom.neighbor_indices = node.adjacent_node_list
	# node1 (n2_old_num)
	atom = dn.dualstruc.atom_list[n2_old_num]
	node = dn.dual_node_list[n2_old_num]
	atom.type = node.type
	atom.neighbor_indices = node.adjacent_node_list
	# node1 (n1_new_num)
	atom = dn.dualstruc.atom_list[n1_new_num]
	node = dn.dual_node_list[n1_new_num]
	atom.type = node.type
	atom.neighbor_indices = node.adjacent_node_list
	# node1 (n2_new_num)
	atom = dn.dualstruc.atom_list[n2_new_num]
	node = dn.dual_node_list[n2_new_num]
	atom.type = node.type
	atom.neighbor_indices = node.adjacent_node_list
	## check the geometrically possible outcome
	ring_quadruple = [dn.dual_node_list[n1_old_num].type
	                  dn.dual_node_list[n2_old_num].type
					  dn.dual_node_list[n1_new_num].type
					  dn.dual_node_list[n2_new_num].type ]
	check_list = Vector{Int64}()
	for num in ring_quadruple
		check = 0
		for num_allowed in dn.r_list
			if num == num_allowed
				check = 1
			end
		end
		push!(check_list,check)
	end
	check_output = true
	for check in check_list
		if check == 0
			check_output = false
		end
	end
	return check_output
end

## find a bond by its nodes
function find_bond_by_nodes(dn::Dual_network,node1,node2)
	for i in 1:length(dn.dual_bond_list)
		bond = dn.dual_bond_list[i]
		if ( (bond.node1.number == node1.number && bond.node2.number == node2.number) ||
			 (bond.node1.number == node2.number && bond.node2.number == node1.number)     )
			return bond
			break
		end
	end
end

## save dual network
function save_dual_network(dn::Dual_network, filename::String)
	output_file = open(filename, "w")
	write(output_file, "# dual network file\n")
	write(output_file, "# num_x num_y\n")
	write(output_file, string(dn.num_x)*" "*string(dn.num_y)*"\n")
	write(output_file, "# bond length:\n")
	write(output_file, string(dn.bond_length)*"\n")
	write(output_file, "# lx ly\n")
	write(output_file, string(dn.lx)*" "*string(dn.ly)*"\n")
	# writing dual nodes
	write(output_file, "# dual nodes:\n")
	write(output_file, "# number of nodes:\n")
	write(output_file, string(length(dn.dual_node_list))*"\n")
	for node in dn.dual_node_list
		write(output_file, string(node.number)*" "*string(node.type)*" "*string(node.pos_x)*" "*string(node.pos_y)*"\n")
	end
	# writing dual bonds
	write(output_file, "# dual bonds:\n")
	write(output_file, "# number of dual bonds:\n")
	write(output_file, string(length(dn.dual_bond_list))*"\n")
	for bond in dn.dual_bond_list
		write(output_file, string(bond.number)*" "*string(bond.node1.number)*" "*string(bond.node2.number)*"\n")
	end
	# writing dual triplets
    write(output_file, "# dual triplets\n")
	write(output_file, "# number of dual triplets: ")
	write(output_file, string(length(dn.dual_triplet_list))*"\n")
	for triplet in dn.dual_triplet_list
		write(output_file, string(triplet.number)*" "*string(triplet.node1.number)*" "*string(triplet.node2.number)*" "*string(triplet.node3.number)*" "*string(triplet.bond1.number)*" "*string(triplet.bond2.number)*" "*string(triplet.bond3.number)*"\n")
	end
	close(output_file)
end

## read dual network
# Nd ... dual nodes; Nb ... dual bonds; Nt ... dual triplets
# nx ... num_x; ny ... num_y
# Nd = nx * ny; Nb = 3 * Nd; Nt = 2 * Nd
#
function read_dual_network(dn::Dual_network, filename::String)
	input_file = open(filename, "r")
	lines = readlines(input_file)
	close(input_file)
	dn.num_x, dn.num_y = readdlm(IOBuffer(lines[3]))
	dn.bond_length = parse(Float64, lines[5])
	dn.lx, dn.ly = readdlm(IOBuffer(lines[7]))
	# writing dual nodes
	dn.dual_node_list = Vector()
	node_list_length = parse(Int, lines[10], base=10)
	for i in 1:node_list_length
		node_index = i + 10
		line_vec = readdlm(IOBuffer(lines[node_index]))
		num = Int(line_vec[1])
		type = Int(line_vec[2])
		pos_x = Float64(line_vec[3])
		pos_y = Float64(line_vec[4])
		push!(dn.dual_node_list, Dual_node(num, type, pos_x, pos_y))
	end
	# writing dual bonds
	dn.dual_bond_list = Vector()
	bond_list_length = parse(Int, lines[node_list_length+13], base=10)
	for ii in 1:bond_list_length
		bond_index = ii + 13 + node_list_length
		line_vec = readdlm(IOBuffer(lines[bond_index]))
		#
		num = Int(line_vec[1])
		node1_num = Int(line_vec[2])
		node2_num = Int(line_vec[3])
		node1 = dn.dual_node_list[node1_num]
		node2 = dn.dual_node_list[node2_num]
		push!(dn.dual_bond_list, Dual_bond(num, node1, node2))
	end
	# writing dual triplets
    dn.dual_triplet_list = Vector()
	st1, st2 = split(lines[node_list_length+bond_list_length+15], ":")
	triplet_list_length = parse(Int, st2, base=10)
	for iii in 1:triplet_list_length
		triplet_index = iii + 15 + node_list_length + bond_list_length
		line_vec = readdlm(IOBuffer(lines[triplet_index]))
		num = Int(line_vec[1])
		node1_num = Int(line_vec[2])
		node2_num = Int(line_vec[3])
		node3_num = Int(line_vec[4])
		bond1 = find_bond_by_nodes(dn,dn.dual_node_list[node1_num],dn.dual_node_list[node2_num])
		bond2 = find_bond_by_nodes(dn,dn.dual_node_list[node2_num],dn.dual_node_list[node3_num])
		bond3 = find_bond_by_nodes(dn,dn.dual_node_list[node3_num],dn.dual_node_list[node1_num])
		#
		triplet = Dual_triplet(num, bond1, bond2, bond3)
		#
		triplet.node1 = dn.dual_node_list[node1_num]
		triplet.node2 = dn.dual_node_list[node2_num]
		triplet.node3 = dn.dual_node_list[node3_num]
		find_triangle_nodes(triplet)
		#
		push!(dn.dual_triplet_list, triplet)
	end
	## get adjacent neighboring nodes
	for i in 1:length(dn.dual_node_list)
		get_adjacent_neigh_nodes(dn.dual_node_list[i], dn.dual_node_list, dn.dual_bond_list)
	end
	##
	build_dual_molstruc(dn)
	## find the neighboring triangles
	for i in 1:length(dn.dual_triplet_list)
		find_neigh_triplets(dn.dual_triplet_list[i],dn.dual_triplet_list)
	end
	## ring calculation
	# calculate neighboring triplets to every node
	for i in 1:length(dn.dual_node_list)
		find_adjacent_neigh_triplets_ordered(dn.dual_node_list[i],dn.dual_triplet_list)
	end
end
