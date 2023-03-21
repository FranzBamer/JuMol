#########################################################################
########### calculate the non-affine displacement field of 2 Boxes
########### input: Box_before, Box_after
#########################################################################

import Core.Array

include("../structure/Structure.jl")

mutable struct Disp_field
	Disp_field() = new()
    deform_1
    deform_2
    disp_field
    aff_disp_field
    non_aff_disp_field
    total_affine_disp1
	total_affine_disp2
	mapped_coor1
	mapped_coor2
end

## initialize the elements of Disp_field
function initialize_disp_field(disp::Disp_field, deform_1, deform_2)
	disp.deform_1 = deform_1
	disp.deform_2 = deform_2
    disp.disp_field = zeros(disp.deform_1.structure.noa,3)
    disp.aff_disp_field = zeros(disp.deform_1.structure.noa,3)
    disp.non_aff_disp_field = zeros(disp.deform_1.structure.noa,3)
    disp.total_affine_disp1 = zeros(disp.deform_1.structure.noa,3)
	disp.total_affine_disp2 = zeros(disp.deform_2.structure.noa,3)
	disp.mapped_coor1 = zeros(disp.deform_1.structure.noa,3)
	disp.mapped_coor2 = zeros(disp.deform_1.structure.noa,3)
end

## Transform to rectangle (total shear deformation)
function transform_to_rectangle_shear_xy(disp::Disp_field)
    tan_gamma = disp.deform_1.structure.box.lxy/disp.deform_1.structure.box.ly
    for i in 1:disp.deform_1.structure.noa
		disp.total_affine_disp[i,1] = tan_gamma * disp.deform_1.structure.atom_list[i].pos[2]
	end
end

function transform_to_rectangle_tension_xy(disp::Disp_field, structure_init)
	for i in 1:disp.deform_1.structure.noa
		disp.total_affine_disp1[i,1] = (disp.deform_1.structure.box.lx - structure_init.box.lx) / structure_init.box.lx * structure_init.atom_list[i].pos[1]
		disp.total_affine_disp1[i,2] = (disp.deform_1.structure.box.ly - structure_init.box.ly) / structure_init.box.ly * structure_init.atom_list[i].pos[2]
		disp.mapped_coor1[i,1] = disp.deform_1.structure.atom_list[i].pos[1] - disp.total_affine_disp1[i,1]
		disp.mapped_coor1[i,2] = disp.deform_1.structure.atom_list[i].pos[2] - disp.total_affine_disp1[i,2]
	end

	for i in 1:disp.deform_2.structure.noa
		disp.total_affine_disp2[i,1] = (disp.deform_2.structure.box.lx - structure_init.box.lx) / structure_init.box.lx * structure_init.atom_list[i].pos[1]
		disp.total_affine_disp2[i,2] = (disp.deform_2.structure.box.ly - structure_init.box.ly) / structure_init.box.ly * structure_init.atom_list[i].pos[2]
		disp.mapped_coor2[i,1] = disp.deform_2.structure.atom_list[i].pos[1] - disp.total_affine_disp2[i,1]
		disp.mapped_coor2[i,2] = disp.deform_2.structure.atom_list[i].pos[2] - disp.total_affine_disp2[i,2]
	end
end

##
function calc_affine_disp_field_simple_shear_xy(disp::Disp_field)
	delta_xy = disp.deform_2.structure.box.lxy - disp.deform_1.structure.box.lxy
	tan_gamma = delta_xy/disp.deform_1.structure.box.ly
	for i in 1:disp.deform_1.structure.noa
	    disp.aff_disp_field[i,1] = tan_gamma * disp.deform_1.structure.atom_list[i].pos[1]
	end
end

##
function calc_affine_disp_field_tension_x(disp::Disp_field,center_x=0.0)
	delta_x = disp.deform_2.structure.box.lx - disp.deform_1.structure.box.lx
	gamma_x = delta_x/disp.deform_1.structure.box.lx
	println(center_x)
	for i in 1:disp.deform_1.structure.noa
		if center_x==0.0
				disp.aff_disp_field[i,1] = gamma_x * disp.deform_1.structure.atom_list[i].pos[1]
		else
	    	disp.aff_disp_field[i,1] = gamma_x * (disp.deform_1.structure.atom_list[i].pos[1]-center_x)
		end
	    disp.aff_disp_field[i,2] = 0.0
	end
end

##
function calc_affine_disp_field_tension_y(disp::Disp_field,center_y=0.0)
	delta_y = disp.deform_2.structure.box.ly - disp.deform_1.structure.box.ly
	gamma_y = delta_y/disp.deform_1.structure.box.ly
	println(center_y)
	for i in 1:disp.deform_1.structure.noa
		if center_y==0.0
				disp.aff_disp_field[i,2] = gamma_y * disp.deform_1.structure.atom_list[i].pos[2]
		else
	    	disp.aff_disp_field[i,2] = gamma_y * (disp.deform_1.structure.atom_list[i].pos[2]-center_y)
		end
	    disp.aff_disp_field[i,1] = 0.0
	end
end

##
function calc_non_aff_disp_field(disp::Disp_field)
	for i in 1:disp.deform_1.structure.noa
		disp.disp_field[i,1] = disp.deform_2.structure.atom_list[i].pos[1] - disp.deform_1.structure.atom_list[i].pos[1]
		disp.disp_field[i,2] = disp.deform_2.structure.atom_list[i].pos[2] - disp.deform_1.structure.atom_list[i].pos[2]
		disp.non_aff_disp_field[i,1] = disp.disp_field[i,1] - disp.aff_disp_field[i,1]
		disp.non_aff_disp_field[i,2] = disp.disp_field[i,2] - disp.aff_disp_field[i,2]
		#disp.disp_field[i,1] = disp.mapped_coor2[i,1] - disp.mapped_coor1[i,1]
		#disp.disp_field[i,1] = disp.mapped_coor2[i,2] - disp.mapped_coor1[i,2]

	end
    #disp.disp_field = get_atom_pos_mat(disp.deform_2.structure)
   	#				- get_atom_pos_mat(disp.deform_1.structure)
    #disp.non_aff_disp_field = disp.disp_field - disp.aff_disp_field
end

##
function erase_pbc_disp(disp::Disp_field)
	lx = disp.deform_1.structure.box.lx
	ly = disp.deform_1.structure.box.ly
	lz = disp.deform_1.structure.box.lz
	for i in 1:disp.deform_1.structure.noa
	    x_dist = abs(disp.disp.non_aff_disp_field[i,1])
		y_dist = abs(disp.disp.non_aff_disp_field[i,2])
		z_dist = abs(disp.disp.non_aff_disp_field[i,3])
	end
	if x_dist > lx * 0.4 | y_dist > ly * 0.4 | z_dist > lz * 0.4
		disp.non_aff_disp_field[i,1] = 0.0
		disp.non_aff_disp_field[i,2] = 0.0
		disp.non_aff_disp_field[i,3] = 0.0
	end
end
