## Structure of the ensemble
######################################
######################################
######################################
#
#
# (c) Franz Bamer, Ivy Wu Mai-2022
######################################



import Core.Array

using DelimitedFiles
using LinearAlgebra

include("Atom.jl")
include("Box.jl")
include("Read.jl")
include("Cell.jl")


mutable struct Structure
	## cutoff radius and skin
	rc::Float64
	rskin::Float64
	## periodic boundary conditions
	pbx::Int
	pby::Int
	pbz::Int
	# linked cell boolean
	linked_cells_bool::Bool
	## load file reader
	reader::Read
	## number of atoms
	noa::Int64
	## atom list
	atom_list::Vector
	## linked cell lists and arrays
	cell_list::Vector
	cell_neighbour_array::Array
	## Box dimensions
	box::Box
	## constructor
	Structure(rc=10.0, rskin=1.0,
			  pbx=0, pby=0, pbz=0,
			  linked_cells_bool=true) = new(rc, rskin,
											pbx, pby, pbz,
											linked_cells_bool)
end

##############################################
#### create empty basic objects for every structure
##############################################
## initialize empty objects of the atomic structure
function initialize_structure_objects(structure::Structure)
	structure.noa = 0
	structure.reader = Read()
	structure.atom_list = Vector()
	structure.box = Box()
	#structure.cell_neighbour_array = Vector()
end

##############################################
#### data loading and reading
##############################################
## load the box from a lammpstrj text file
function read_lammpstrj(structure::Structure,filename,timestep_start=1, triclinic=false)
	initialize_structure_objects(structure)
	get_lines(structure.reader,filename)
	read_lammpstrj(structure.reader, structure, structure.reader.line_start_list[timestep_start], triclinic)
	set_box_basis_vectors(structure.box)
end

## load the box from a lammps type box input file
function read_lammps_box(structure::Structure,filename)
	initialize_structure_objects(structure)
	read_box_lammps(structure)
end

## load an atom by hand
function add_atom_by_hand(structure::Structure,number,type,x,y,z,vx,vy,vz,group=0)
	push!(structure.atom_list, Atom(number,type))
	structure.atom_list[structure.noa+1].pos = [x,y,z]
	structure.atom_list[structure.noa+1].vel = [vx,vy,vz]
	structure.atom_list[structure.noa+1].acc = [0.0,0.0,0.0]
	structure.atom_list[structure.noa+1].group = group
	structure.noa += 1
end

## load the box dimensions by hand
function create_box_by_hand(structure::Structure,lx,ly,lz,lxy,lyz,lxz)
	structure.box.lx = lx
	structure.box.ly = ly
	structure.box.lz = lz
	structure.box.lxy = lxy
	structure.box.lyz = lyz
	structure.box.lxz = lxz
	set_box_basis_vectors(structure.box)
end

## get the positions in a matrix 3 x noa
function get_atom_pos_mat(structure::Structure)
	atom_pos_mat = zeros(structure.noa,3)
	for i in 1:structure.noa
		atom_pos_mat[i,:] = structure.atom_list[i].pos[1:3]
	end
	return atom_pos_mat
end

## get global position vector noa*3 x 1
function get_global_atom_pos_vec(structure::Structure)
	atom_pos_vec = zeros(structure.noa*3)
	cntr = 1
	for i in 1:structure.noa
		for ii in 1:3
			atom_pos_vec[cntr] = structure.atom_list[i].pos[ii]
			cntr += 1
		end
	end
	return atom_pos_vec
end

## set atomic positions according to a global position vector
function set_global_pos_vec(structure::Structure, pos_vec)
	cntr = 1
	for i in 1:structure.noa
		for ii in 1:3
			structure.atom_list[i].pos[ii] = pos_vec[cntr]
			cntr += 1
		end
	end
end

## get global force vector noa*3 x 1
function get_global_force_vec(structure::Structure)
	force_vec = zeros(structure.noa*3)
	cntr = 1
	for i in 1:structure.noa
		for ii in 1:3
			force_vec[cntr] = structure.atom_list[i].force[ii]
			cntr += 1
		end
	end
	return force_vec
end

## set atomic forces according to a global force vector
function set_global_force_vec(structure::Structure, force_vec)
	cntr = 1
	for i in 1:structure.noa
		for ii in 1:3
			structure.atom_list[i].force[ii] = force_vec[cntr]
			cntr += 1
		end
	end
end

##############################################
#### distance and neighbor list calculations
##############################################
## initialize structure (must be run once)
function initialize_structure(structure::Structure, include_dist_updates=true)
	## create matrices for distances for every atom (4,noa-num_atom)
	for i in 1:structure.noa
		structure.atom_list[i].distances = zeros(4,structure.noa-i)
	end
	## initialize cell lists
	update_cell(structure)
	linked_cell_list(structure)
	construct_neigh_cell_list(structure)
	## full distance update
	if include_dist_updates
		update_distances(structure)
		update_neighbor_lists(structure)
	end
end

## put atoms back to box if they leave the box
function put_atoms_back_to_box(structure::Structure)
	for i in 1:structure.noa
		## to do (transformation into the triclinic coordinate system needed!!!!!!)
		if structure.atom_list[i].pos[1] > (structure.box.h1[1] +
			                                structure.box.h2[1]/structure.box.h2[2]* structure.atom_list[i].pos[2] +
											structure.box.h3[1]/structure.box.h3[3]*structure.atom_list[i].pos[3]   )
			structure.atom_list[i].pos -= structure.box.h1
		elseif structure.atom_list[i].pos[1] < 0# +
			                                    #structure.box.h2[1]/structure.box.h2[2]* structure.atom_list[i].pos[2] +
			                                    #structure.box.h3[1]/structure.box.h3[3]*structure.atom_list[i].pos[3]   ) #FIXME check if this is correct!!!!!??????
			structure.atom_list[i].pos += structure.box.h1
		end
		if structure.atom_list[i].pos[2] > structure.box.h2[2]
			structure.atom_list[i].pos -= structure.box.h2
		elseif structure.atom_list[i].pos[2] < 0
			structure.atom_list[i].pos += structure.box.h2
		end
		if structure.atom_list[i].pos[3] > structure.box.h3[3]
			structure.atom_list[i].pos -= structure.box.h3
		elseif structure.atom_list[i].pos[3] < 0
			structure.atom_list[i].pos += structure.box.h3
		end
	end
end

## inlcude distance update functions
include("Distance.jl")
