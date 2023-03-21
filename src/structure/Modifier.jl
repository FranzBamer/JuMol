	# Atomic structure
######################################
######################################
######################################


include("Atom.jl")


mutable struct Modifier
	structure
	Modifier(structure) = new(structure)
end


## translate a sample over periodic boundary conditions
function translate(mod::Modifier, translate_vector)
	for atom in mod.structure.atom_list
		atom.pos[1:2] += translate_vector
	end
end

## rotate circular sample by a certain angle counter-clockwise
function rotate(mod::Modifier, theta::Float64)
	rot_mat = [cosd(theta) -sind(theta)
	           sind(theta)  cosd(theta)]
	for atom in mod.structure.atom_list
		atom.pos[1:2] = rot_mat*atom.pos[1:2]
	end
end

## cut circle from box
function cut_circle_from_box(mod::Modifier,center_x,center_y,radius)
    new_atom_list = Vector()
    cntr = 1
    for i in 1:mod.structure.noa
		type = mod.structure.atom_list[i].type
        pos_x = mod.structure.atom_list[i].pos[1]
        pos_y = mod.structure.atom_list[i].pos[2]
        pos_z = mod.structure.atom_list[i].pos[3]
        dist = sqrt( (pos_x-center_x)^2 + (pos_y-center_y)^2 )
        if dist < radius
			atom = Atom(cntr,type)
			atom.pos = copy(mod.structure.atom_list[i].pos)
			atom.vel = copy(mod.structure.atom_list[i].vel)
            atom.group = mod.structure.atom_list[i].group
			push!(new_atom_list,atom)
			cntr += 1
        end
    end
    mod.structure.atom_list = new_atom_list
    mod.structure.noa = length(new_atom_list)
    mod.structure.box.lx = 2*radius
    mod.structure.box.ly = 2*radius
	mod.structure.box.lxy = 0.0
	mod.structure.box.lxz = 0.0
	mod.structure.box.lyz = 0.0
	set_box_basis_vectors(mod.structure.box)
end

## define group out of a circle
function define_group_out_of_circle(mod::Modifier,center_x,center_y,radius)
    for i in 1:mod.structure.noa
        pos_x = mod.structure.atom_list[i].pos[1]
        pos_y = mod.structure.atom_list[i].pos[2]
        pos_z = mod.structure.atom_list[i].pos[3]
        dist = sqrt( (pos_x-center_x)^2 + (pos_y-center_y)^2 )
        if dist > radius
            mod.structure.atom_list[i].group = 1
        end
    end
end

## rotation region
function rotate_circle(mod::Modifier,center_x,center_y,center_z,axis,degree)
	if axis=='x'
		mat_rot = [[1.0,0,0] [0,cosd(degree),sind(degree)] [0,-sind(degree),cosd(degree)]]
	end
	if axis=='y'
		mat_rot = [[cosd(degree),0,-sind(degree)] [0,1.0,0] [sind(degree),0,cosd(degree)]]
	end
	if axis=='z'
		mat_rot = [[cosd(degree),sind(degree),0] [-sind(degree),cosd(degree),0] [0,0,1.0]]
	end

	center_vec = [center_x,center_y,center_z]
	for i in 1:mod.structure.noa
		mod.structure.atom_list[i].pos = mat_rot * (mod.structure.atom_list[i].pos - center_vec) + center_vec
	end

end


#### functions to check phenomena
## shuffle the atom indices
function shuffle_atom_indices(mod::Modifier)
	# shuffle
	mod.structure.atom_list = shuffle(mod.structure.atom_list)
	# create new list
	new_atom_list = Vector()
    cntr = 1
    for atom in mod.structure.atom_list
		type = atom.type
        pos_x = atom.pos[1]
        pos_y = atom.pos[2]
        pos_z = atom.pos[3]
		atom_new = Atom(cntr, type)
		atom_new.pos = copy(atom.pos)
		atom_new.vel = copy(atom.vel)
        atom_new.group = atom.group
		push!(new_atom_list, atom_new)
		cntr += 1
    end
    mod.structure.atom_list = new_atom_list
end
