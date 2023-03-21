## Deformations
######################################
######################################
######################################
## methods to deform the box and atoms within the box

mutable struct Aff_deform
    structure
    Aff_deform(structure) = new(structure)
	delta_x::Float64
	delta_y::Float64
	delta_z::Float64
	delta_xy::Float64
	delta_xz::Float64
	delta_yz::Float64
end

function set_affine_deform_vol(deform::Aff_deform,
	                           delta_x::Float64,
							   delta_y::Float64,
							   delta_z::Float64)
	deform.delta_x = delta_x
	deform.delta_y = delta_y
	deform.delta_z = delta_z
	## affine atom deformation
	for i in 1:deform.structure.noa
		deform.structure.atom_list[i].pos[1] += delta_x/deform.structure.box.lx * deform.structure.atom_list[i].pos[1]
		deform.structure.atom_list[i].pos[2] += delta_y/deform.structure.box.ly * deform.structure.atom_list[i].pos[2]
		deform.structure.atom_list[i].pos[3] += delta_z/deform.structure.box.lz * deform.structure.atom_list[i].pos[3]
	end
	## deform box
	deform.structure.box.lx += delta_x
	deform.structure.box.ly += delta_y
	deform.structure.box.lz += delta_z
	set_box_basis_vectors(deform.structure.box)
end

function set_affine_deform_shear(deform::Aff_deform,delta_xy,delta_xz,delta_yz)
	deform.delta_xy = delta_xy
	deform.delta_xz = delta_xz
	deform.delta_yz = delta_yz
	## affine deformation
	for i in 1:deform.structure.noa
		deform.structure.atom_list[i].pos[1] += delta_xy/deform.structure.box.ly * deform.structure.atom_list[i].pos[2]
		deform.structure.atom_list[i].pos[1] += delta_xz/deform.structure.box.lz * deform.structure.atom_list[i].pos[3]
		deform.structure.atom_list[i].pos[2] += delta_yz/deform.structure.box.lz * deform.structure.atom_list[i].pos[3]
	end
	## deform box
	deform.structure.box.lxy += delta_xy
	deform.structure.box.lxz += delta_xz
	deform.structure.box.lyz += delta_yz
	set_box_basis_vectors(deform.structure.box)
end

function set_true_shear_to_atom_group_2d(deform::Aff_deform,center_x,center_y,delta_x)
	deform.delta_x = delta_x
    delta_y = deform.structure.box.lx * deform.structure.box.ly / (deform.structure.box.lx+delta_x) - deform.structure.box.ly
	deform.delta_y = delta_y
	for i in 1:deform.structure.noa
        deform.structure.atom_list[i].pos[1] = (1 + delta_x/deform.structure.box.lx) * (deform.structure.atom_list[i].pos[1] - center_x) + center_x
        deform.structure.atom_list[i].pos[2] = (1 + delta_y/deform.structure.box.ly) * (deform.structure.atom_list[i].pos[2] - center_y) + center_y
    end
    ## deform box
	#deform.structure.box.ly = deform.structure.box.lx*deform.structure.box.ly/(deform.structure.box.lx+delta_x)
	deform.structure.box.lx += delta_x
    deform.structure.box.ly += delta_y
    set_box_basis_vectors(deform.structure.box)
end

function set_tension_to_atom_group_2d(deform::Aff_deform,center_x,center_y,delta_x,delta_y)
	deform.delta_x = delta_x
	deform.delta_y = delta_y
	for i in 1:deform.structure.noa
        deform.structure.atom_list[i].pos[1] = (1 + delta_x/deform.structure.box.lx) * (deform.structure.atom_list[i].pos[1] - center_x) + center_x + delta_x/2
        deform.structure.atom_list[i].pos[2] = (1 + delta_y/deform.structure.box.ly) * (deform.structure.atom_list[i].pos[2] - center_y) + center_y + delta_y/2
    end
    ## deform box
	deform.structure.box.lx += delta_x
    deform.structure.box.ly += delta_y
    set_box_basis_vectors(deform.structure.box)
end
