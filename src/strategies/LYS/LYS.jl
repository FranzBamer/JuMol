## struct for LYS method
######################################
######################################
######################################
#
# Input: LI (mutable struct), CI (mutable struct), SI (mutable struct), r_free::Float64
#        #
#        LI must contain:
#        (r_free_list::Vector), num_x_grid::Int64, num_y_grid::Int64,
#        delta_theta::Float64, output_path::String, r_total::Float64
#        #
#        CI must contain:
#        (num_samples::Int64), include_lmp::Bool, save_video::Bool,
#        lmp_executable::String, rc::Float64, skin::Float64
#        #
#        SI must contain:
#        (num_x::Int64), (num_y::Int64), (sig::Int64),
#        input_path::String, filename::String
#        #
#        PI should contain (if plotting is performed):
#        hfig::Int64, bfig::Int64, box_bound_col::String,
#        box_bound_lw::Float64, at_size1::Float64, at_size2::Float64,
#        at_col1::String, at_col2::String, at_lw::Float64,
#        factor_vec_x::Vector{Float64}, factor_vec_y::Vector{Float64}
#
# (c) Franz Bamer, Ivy Wu Mai-2022
######################################

# this file must be generalized to local probing in general
# local potential energy, local free volume, local ring statistics,
# local ring neighborhood statistics, local fabric tensors,

## include struct that analyzes the local structure
include("Local_structure.jl")

## LYS method
mutable struct LYS
	LI #(struct defined by input)
	CI #(struct defined by input)
	SI #(struct defined by input)
	PI #(struct defined by input)
	r_free::Float64
	num_sample::Int64
	deform_type::Int64 # 1 for true shear; 2 for simple tension
    LYS(LI, CI, SI, PI, r_free, num_sample, deform_type=1) = new(LI, CI, SI, PI, r_free, num_sample, deform_type)
	## calculation properties (must be defined)
	molstruc
	## output properties
	file_identification::String
	## scanning parameters
	x_coord_list::Vector
	y_coord_list::Vector
	theta_list::Vector
end


## load the sample
function load_sample(lys::LYS)
	file_load = lys.SI.input_path*lys.SI.filename
	println("... importing sample ...")
	println(file_load)
	## load entire ensemble
	lys.molstruc = Jumol.Structure()
	lys.molstruc.rc = lys.CI.rc # pseudo cutoff
	lys.molstruc.rskin = lys.CI.rskin # pseudo cutoff skin
	lys.molstruc.pbx=1; lys.molstruc.pby=1; # periodic boundary conditions
	initialize_structure_objects(lys.molstruc)
	## read lammpstrj
	read_lammpstrj(lys.molstruc, file_load)
	initialize_structure(lys.molstruc)
	## create file identification
	lys.file_identification = "sig"*string(lys.SI.sig)*"_num"*string(lys.num_sample)*"_r"*string(Int(lys.r_free))
end


## plot the sample
function plot_sample(lys::LYS)
	println("... plotting sample ", lys.num_sample, " ...")
    Plotter = Vis2d(lys.molstruc)
    Plotter.hfig = lys.PI.hfig
    Plotter.bfig = lys.PI.bfig
    fig = plot_box(Plotter, lys.PI.box_bound_col, lys.PI.box_bound_lw)
	lx0 = lys.molstruc.box.lx; ly0 = lys.molstruc.box.ly
    plot_atomic_structure_binary(Plotter, lys.PI.at_size1, lys.PI.at_size2,
	                             lys.PI.factor_vec_x*lx0,
			 		             lys.PI.factor_vec_x*ly0,
								 lys.PI.at_col1, lys.PI.at_col2, lys.PI.at_lw)
    display(fig)
	savefig(fig,lys. LI.output_path*"sample_"*lys.file_identification*".pdf")
end

## run LYS function
function run_lys(lys::LYS, x_coord_start=0.0, y_coord_start=0.0)
	println("... running LYS scan ...")
	println("r_free: ", lys.r_free)
	## generate the grid and prepare lists
	lys.x_coord_list = Vector(); lys.y_coord_list = Vector(); lys.theta_list = Vector();
	x_coord = x_coord_start; y_coord = y_coord_start; theta = 0.0;
	dx = lys.molstruc.box.lx/lys.LI.num_x_grid; dy = lys.molstruc.box.ly/lys.LI.num_y_grid
	num_angle = Int(180.0/lys.LI.delta_theta)
	theta = -90.0
	for i in 1:lys.LI.num_x_grid
		push!(lys.x_coord_list, x_coord)
		x_coord += dx
	end
	for i in 1:lys.LI.num_y_grid
		push!(lys.y_coord_list, y_coord)
		y_coord += dy
	end
	for i in 1:num_angle
		push!(lys.theta_list, theta)
		theta += lys.LI.delta_theta
	end
	## create stress_strain_output_file
	stress_strain_file = open(lys.LI.output_path*"stress_strain_lys_output_"*lys.file_identification*".dat", "w")
	## go through the grid
	cntr_x = 1
	for x_coord in lys.x_coord_list
		cntr_y = 1
		for y_coord in lys.y_coord_list
			cntr_angle = 1
			for theta in lys.theta_list
				println("cntr_x: ", cntr_x, ", cntr_y: ", cntr_y)
				## create new molstruc for the local description
				loc_molstruc = Jumol.Structure()
				##
				loc_molstruc.rc = lys.CI.rc; loc_molstruc.rskin = lys.CI.rskin; loc_molstruc.pbx = 0; loc_molstruc.pby = 0
				Jumol.initialize_structure_objects(loc_molstruc)
				Jumol.create_box_by_hand(loc_molstruc, lys.molstruc.box.lx, lys.molstruc.box.ly, lys.molstruc.box.lz,
				                                       lys.molstruc.box.lxy, lys.molstruc.box.lyz, lys.molstruc.box.lxz)
				for atom in lys.molstruc.atom_list
					Jumol.add_atom_by_hand(loc_molstruc, atom.number, atom.type,
					                       atom.pos[1], atom.pos[2], 0.0,
										   atom.vel[1], atom.vel[2], 0.0)
				end
				# move box from the original center coordinate to the grid point
				# so that the grid is the center considering periodic boundary
				# conditions
				mol_modifier = Jumol.Modifier(loc_molstruc)
				local_point = [x_coord
				               y_coord]
				center_point = [lys.molstruc.box.lx*0.5
				                lys.molstruc.box.ly*0.5]
				translate_vector = center_point - local_point
				Jumol.translate(mol_modifier, translate_vector)
				Jumol.put_atoms_back_to_box(loc_molstruc)
				# cutout a circle around the new center with the total radius
				Jumol.cut_circle_from_box(mol_modifier, loc_molstruc.box.lx*0.5,loc_molstruc.box.ly*0.5, lys.LI.r_total)
				# translate center of sample into the origin of the
				# coordinate system
				Jumol.translate(mol_modifier, center_point*(-1.0))
				# rotate the sample by theta
				Jumol.rotate(mol_modifier, theta)
				# translate origin into left bottom box corner
				sample_center = [loc_molstruc.box.lx*0.5
				                 loc_molstruc.box.ly*0.5]
				Jumol.translate(mol_modifier, sample_center)
				# define the atom group outside r_free to 1 so that the atoms can
				# be fixed in this area during the calculation
				Jumol.define_group_out_of_circle(mol_modifier, sample_center[1], sample_center[2], lys.r_free)
				Jumol.initialize_structure(loc_molstruc)
				## load local object for the mechanical investigation
				Local_point = Local_structure(loc_molstruc, lys.r_free, cntr_x, cntr_y, x_coord, y_coord, theta, lys.LI.output_path, lys.CI.save_video)
				if lys.CI.include_lmp
                    if lys.deform_type ==1
    					run_true_shear_deformation_lmp(Local_point, lys.CI.lmp_executable, stress_strain_file, lys.file_identification)
                    else
                        run_simple_tension_deformation_lmp(Local_point, lys.CI.lmp_executable, stress_strain_file, lys.file_identification)
                    end
				else
                    if lys.deform_type ==1
    					run_true_shear_deformation(Local_point, stress_strain_file, cntr_angle)
                    else
                        run_simple_tension_deformation(Local_point, stress_strain_file, cntr_angle)
                    end
				end
				cntr_angle += 1
			end
			cntr_y += 1
		end
		cntr_x += 1
	end
	## close stress_strain_output_file
	close(stress_strain_file)
end
