## struct for the LYS method
######################################
######################################
######################################
#
#
# (c) Franz Bamer, Ivy Wu Mai-2022
######################################

include("Local_structure.jl")


mutable struct LYS
	output_path::String
	num_x_grid::Int64
	num_y_grid::Int64
	dtheta::Float64
	r_total::Float64
	r_free::Float64
	include_lmp::Bool
	save_local_video::Bool
	LYS(output_path,
		num_x_grid, num_y_grid, dtheta,
		                r_total, r_free,
						include_lmp=false,
						save_local_video=false) = new(output_path,
											    num_x_grid, num_y_grid, dtheta,
											    r_total, r_free,
												include_lmp,
												save_local_video)
	## main file properties
	input_path::String
	filename::String
	lmp_executable::String
	molstruc
	## scanning parameters
	x_coord_list
	y_coord_list
	theta_list

end

## load the sample
function load_sample(lys::LYS, input_path, filename)
	lys.input_path = input_path
	lys.filename = filename
	file_load = lys.input_path*lys.filename
	println("... importing molecular file ...")
	println(file_load)
	## load entire ensemble
	lys.molstruc = Jumol.Structure()
	lys.molstruc.rc=2.0 # pseudo cutoff
	lys.molstruc.rskin=0.4 # pseudo cutoff skin
	lys.molstruc.pbx=1; lys.molstruc.pby=1;
	initialize_structure_objects(lys.molstruc)
	## read lammpstrj
	read_lammpstrj(lys.molstruc, file_load)
	initialize_structure(lys.molstruc)
end


## load the sample from Spencer #FIXME #erase later
function load_sample_spencer(lys::LYS, input_path, filename)
	lys.input_path = input_path
	lys.filename = filename
	file_load = lys.input_path*lys.filename
	println("... importing molecular file ...")
	println(file_load)
	## load entire ensemble
	lys.molstruc = Jumol.Structure()
	lys.molstruc.rc=2.0 # pseudo cutoff
	lys.molstruc.rskin=0.4 # pseudo cutoff skin
	lys.molstruc.pbx=1; lys.molstruc.pby=1;
	initialize_structure_objects(lys.molstruc)
	## read lammpstrj
	get_lines(lys.molstruc.reader,filename)
	read_lammpstrj_spencer(lys.molstruc.reader, lys.molstruc, lys.molstruc.reader.line_start_list[end])
	set_box_basis_vectors(lys.molstruc.box)
	initialize_structure(lys.molstruc)
end

## plot a sample
function plot_sample(lys::LYS,
	                 size1=5.5, size2=3.5, lw=0.25,
	                 hfig=1000, bfig=1000,
					 col1="blue", col2="red")
    Plotter = Vis2d(lys.molstruc)
    Plotter.hfig = hfig
    Plotter.bfig = bfig
    fig = plot_box(Plotter,"gray",1.0)
	lx0 = lys.molstruc.box.lx; ly0 = lys.molstruc.box.ly
    plot_atomic_structure_binary(Plotter, size1, size2,
	                             [lx0*(-0.05), lx0*(1.05), lx0*(1.05), lx0*(-0.05)],
			 		             [ly0*(-0.05), ly0*(-0.05), ly0*(1.05), ly0*(1.05)],
								 col1, col2,
								 lw)
    display(fig)
	savefig(fig,"sample.pdf")
end

## run LYS function
function run_lys(lys::LYS, x_coord_start=0.0, y_coord_start=0.0)
	## generate the grid
	lys.x_coord_list = Vector(); lys.y_coord_list = Vector(); lys.theta_list = Vector();
	x_coord = x_coord_start; y_coord = y_coord_start; theta = 0.0;
	dx = lys.molstruc.box.lx/lys.num_x_grid; dy = lys.molstruc.box.ly/lys.num_y_grid
	#num_angle = Int(180.0/lys.dtheta) #FIXME
	num_angle = 10
	theta = 0.0
	for i in 1:lys.num_x_grid
		push!(lys.x_coord_list, x_coord)
		x_coord += dx
	end
	for i in 1:lys.num_y_grid
		push!(lys.y_coord_list, y_coord)
		y_coord += dy
	end
	for i in 1:num_angle
		push!(lys.theta_list, theta)
		theta += lys.dtheta
	end
	## create stress_strain_output_file
	stress_strain_file = open(lys.output_path*"stress_strain_lys_output.dat", "w")
	## go through the grid
	cntr_x = 1
	for x_coord in lys.x_coord_list
		cntr_y = 1
		for y_coord in lys.y_coord_list
			cntr_angle = 1
			for theta in lys.theta_list
				## create new molstruc for the local description
				loc_molstruc = Jumol.Structure()
				##FIXME check if the right cutoff is defined!!!!!!
				loc_molstruc.rc = 2.0; loc_molstruc.rskin = 0.2; loc_molstruc.pbx = 0; loc_molstruc.pby = 0
				Jumol.initialize_structure_objects(loc_molstruc)
				Jumol.create_box_by_hand(loc_molstruc, lys.molstruc.box.lx, lys.molstruc.box.ly, lys.molstruc.box.lz,
				                                       lys.molstruc.box.lxy, lys.molstruc.box.lyz, lys.molstruc.box.lxz)
				for atom in lys.molstruc.atom_list
					Jumol.add_atom_by_hand(loc_molstruc, atom.number, atom.type,
					                       atom.pos[1], atom.pos[2], 0.0,
										   atom.vel[1], atom.vel[2], 0.0)
				end
				#
				## modify the loc_molstruc_object
				mol_modifier = Jumol.Modifier(loc_molstruc)
				Jumol.shuffle_atom_indices(mol_modifier) #FIXME erase again later
				println("check shuffling")
				println("position atom 1: ")
				println(loc_molstruc.atom_list[1].pos)
				println("position atom 2: ")
				println(loc_molstruc.atom_list[2].pos)
				local_point = [x_coord
				               y_coord]
				center_point = [lys.molstruc.box.lx*0.5
				                lys.molstruc.box.ly*0.5]
				translate_vector = center_point - local_point
				Jumol.translate(mol_modifier, translate_vector)
				Jumol.put_atoms_back_to_box(loc_molstruc)
				Jumol.cut_circle_from_box(mol_modifier, loc_molstruc.box.lx*0.5,loc_molstruc.box.ly*0.5, lys.r_total)
				# translate center of sample into the origin
				Jumol.translate(mol_modifier, center_point*(-1.0))
				# rotate by theta
				Jumol.rotate(mol_modifier, theta)
				# translate origin into left box corner
				sample_center = [loc_molstruc.box.lx*0.5
				                 loc_molstruc.box.ly*0.5]
				Jumol.translate(mol_modifier, sample_center)
				Jumol.define_group_out_of_circle(mol_modifier, sample_center[1], sample_center[2], lys.r_free)
				Jumol.initialize_structure(loc_molstruc)
				## load local object
				Local_point = Local_structure(loc_molstruc, lys.r_free, cntr_x, cntr_y, x_coord, y_coord, theta, lys.output_path, lys.save_local_video)
				if lys.include_lmp
					run_true_shear_deformation_lmp(Local_point, lys.lmp_executable, stress_strain_file)
				else
					run_true_shear_deformation(Local_point, stress_strain_file, cntr_angle)
				end
				## write the lammpstrj of the local structure #FIXME
				open_file(loc_molstruc.reader, "sample_shuffle_"*string(cntr_angle)*".lammpstr", "w") #FIXME
				write_box(loc_molstruc.reader,loc_molstruc,1) #FIXME
				close_file(loc_molstruc.reader) #FIXME
				cntr_angle += 1
			end
			cntr_y += 1
		end
		cntr_x += 1
	end
	## close stress_strain_output_file
	close(stress_strain_file)
end
