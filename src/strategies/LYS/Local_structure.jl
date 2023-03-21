## struct for the local response
######################################
######################################
######################################
#
#
# (c) Franz Bamer, Ivy Wu Mai-2022
######################################





mutable struct Local_structure
	molstruc
	r_free::Float64
	cntr_x::Int64
	cntr_y::Int64
	x_coord::Float64
	y_coord::Float64
	theta::Float64
	output_path::String
	save_def_vid::Bool
	Local_structure(molstruc, r_free,
	                cntr_x, cntr_y,
	                x_coord, y_coord,
					theta,
					output_path,
					save_def_vid) = new(molstruc, r_free,
					             cntr_x, cntr_y,
				                 x_coord, y_coord,
				                 theta,
								 output_path,
								 save_def_vid)
	lmp_executable::String
end



## run true shear AQS simulation
function run_true_shear_deformation(ls::Local_structure, stress_strain_output_file)
	# get the number of atoms in r_free
	cntr_atom_free = 0
	for atom in ls.molstruc.atom_list
		if atom.group == 0
			cntr_atom_free += 1
		end
	end
	#
	println("########")
	println("Local point scan:")
	println("num_x: ", ls.cntr_x, ", num_y: ", ls.cntr_y)
	println("coord_x: ", ls.x_coord, ", coord_y: ", ls.y_coord)
	println("theta: ", ls.theta)
	println("number of total atoms: ", length(ls.molstruc.atom_list))
	println("number of free atoms: ", cntr_atom_free)
	println("##")
	# video output filename
	video_filename = ls.output_path*"videos/deform_"*string(ls.cntr_x)*"_"*string(ls.cntr_y)*"_"*string(ls.theta)
	# stress strain output file
	#
	deformation_steps = 1000
	#
	molcalc = Jumol.Calc(ls.molstruc)
	Jumol.initialize_potential(molcalc, 13)
	#
	eps = 1.0e-10
	max_num_steps = 100000
	Minimizer = Jumol.Min(ls.molstruc, molcalc)
	Minimizer.alpha_min = 1.0e0
	Minimizer.acc_factor = 1.001
	Jumol.run_cg(Minimizer, max_num_steps, eps)
	#
	# get the initial stress
	Jumol.calc_stress_tensor(molcalc)
	println("stress tensor:")
	pretty_table(molcalc.stress_tensor, noheader = true,
		         crop = :horizontal, formatters = ft_round(8))
	stress_before = molcalc.stress_tensor
	write(stress_strain_output_file, string(stress_before[1,1])*" "*string(stress_before[2,2])*" "*string(stress_before[1,2])*" ")
	#
	deform = Jumol.Aff_deform(ls.molstruc)
	lx0 = ls.molstruc.box.lx
	ly0 = ls.molstruc.box.ly
	delta_x = 0.1
	delta_stress = 1.0
	cntr_step = 1
	while delta_stress > 0.0
		#
		# step number
		println("####")
		println("number of AQS step: ", cntr_step)
		# box deformation in transversal direction
		delta_y = -(ls.molstruc.box.ly*delta_x)/(ls.molstruc.box.lx+delta_x)
		# deform box
		Jumol.set_affine_deform_vol(deform, delta_x, delta_y, 0.0)
		# minimize potential energy
		Minimizer.alpha_min = 1.0e0
		Minimizer.acc_factor = 1.001
		Jumol.run_cg(Minimizer, max_num_steps, eps)
		# plot sample
		plot_sample(ls, ls.molstruc,
		            [lx0*(-0.05), lx0*(1.50), lx0*(1.50), lx0*(-0.05)],
		            [ly0*(-0.05), ly0*(-0.05), ly0*(1.05), ly0*(1.05)],
					7.0, 12.0)
		# get the stress
		Jumol.calc_stress_tensor(molcalc)
		println("stress tensor:")
	    pretty_table(molcalc.stress_tensor, noheader = true,
	        crop = :horizontal, formatters = ft_round(8))
		delta_stress = molcalc.stress_tensor[1,1] - stress_before[1,1]
		# update stress tensor
		stress_before = molcalc.stress_tensor
		#
		if delta_stress < 0
			write(stress_strain_output_file, string(cntr_step)*" "*string(stress_before[1,1])*" "*string(stress_before[2,2])*" "*string(stress_before[1,2])*"\n")
		end
		#println(delta_stress)
		#readline()
		cntr_step += 1
	end
	#if ls.save_def_vid
	#	gif(anim, video_filename*".gif", fps = 10)
	#end
end

## true shear deformation with lammps
function run_true_shear_deformation_lmp(ls::Local_structure, lmp_executable::String, stress_strain_output_file, loc_file_identification)
	# get the number of atoms in r_free
	cntr_atom_free = 0
	for atom in ls.molstruc.atom_list
		if atom.group == 0
			cntr_atom_free += 1
		end
	end
	#
	println("########")
	println("Local point scan:")
	println("num_x: ", ls.cntr_x, ", num_y: ", ls.cntr_y)
	println("coord_x: ", ls.x_coord, ", coord_y: ", ls.y_coord)
	println("theta: ", ls.theta)
	println("number of total atoms: ", length(ls.molstruc.atom_list))
	println("number of free atoms: ", cntr_atom_free)
	println("##")
	loc_sample_identification = "_cntrx"*string(ls.cntr_x)*"_cntry"*string(ls.cntr_y)*"_theta"*string(ls.theta)
	# write lammps input files
	Jumol.open_file(ls.molstruc.reader, "loc_sample_group_0.twodsilica", "w")
	Jumol.write_lmp_box_group(ls.molstruc.reader, ls.molstruc, 1, 0)
	Jumol.close_file(ls.molstruc.reader)
	Jumol.open_file(ls.molstruc.reader, "loc_sample_group_1.twodsilica", "w")
	Jumol.write_lmp_box_group(ls.molstruc.reader, ls.molstruc, 1, 1)
	Jumol.close_file(ls.molstruc.reader)
	# run the deformation with lammps
	run(`$lmp_executable -in deform.in`)
	println("lammps calculation finished")
	# read lammps stress output file
	strain_xx_vec, strain_yy_vec, stress_xx_vec, stress_yy_vec, stress_xy_vec, step_vec = read_lmp_stress_file(ls)
	# stress strain output file
	write(stress_strain_output_file, "#####\n")
	write(stress_strain_output_file, "number of atoms:\n")
	write(stress_strain_output_file, string(length(ls.molstruc.atom_list))*"\n")
	write(stress_strain_output_file, "number of atoms free:\n")
	write(stress_strain_output_file, string(cntr_atom_free)*"\n")
	write(stress_strain_output_file, "num_x_grid, num_y_grid, theta\n")
	write(stress_strain_output_file, string(ls.cntr_x)*" "*string(ls.cntr_y)*" "*string(ls.theta)*"\n")
	write(stress_strain_output_file, "pos_x, pos_y\n")
	write(stress_strain_output_file, string(ls.x_coord)*" "*string(ls.y_coord)*"\n")
	write(stress_strain_output_file, "# init_sigma_xx, init_sigma_yy, init_sigma_xy, sigma_xx, sigma_yy, sigma_xy\n")
	write(stress_strain_output_file, string(stress_xx_vec[1])*" "*string(stress_yy_vec[1])*" "*string(stress_xy_vec[1])*" "*
	                                 string(stress_xx_vec[end-1])*" "*string(stress_yy_vec[end-1])*" "*string(stress_xy_vec[end-1])*"\n")
	# video output filename
	video_output_filename = ls.output_path*"locdef_"*loc_file_identification*loc_sample_identification
	if ls.save_def_vid
		save_deformation_vid(ls, video_output_filename)
	end
end

## run simple tension AQS simulation
function run_simple_tension_deformation(ls::Local_structure, stress_strain_output_file)
	# get the number of atoms in r_free
	cntr_atom_free = 0
	for atom in ls.molstruc.atom_list
		if atom.group == 0
			cntr_atom_free += 1
		end
	end
	#
	println("########")
	println("Local point scan:")
	println("num_x: ", ls.cntr_x, ", num_y: ", ls.cntr_y)
	println("coord_x: ", ls.x_coord, ", coord_y: ", ls.y_coord)
	println("theta: ", ls.theta)
	println("number of total atoms: ", length(ls.molstruc.atom_list))
	println("number of free atoms: ", cntr_atom_free)
	println("##")
	# video output filename
	video_filename = ls.output_path*"videos/deform_"*string(ls.cntr_x)*"_"*string(ls.cntr_y)*"_"*string(ls.theta)
	# stress strain output file
	#
	deformation_steps = 1000
	#
	molcalc = Jumol.Calc(ls.molstruc)
	Jumol.initialize_potential(molcalc, 13)
	#
	eps = 1.0e-10
	max_num_steps = 100000
	Minimizer = Jumol.Min(ls.molstruc, molcalc)
	Minimizer.alpha_min = 1.0e0
	Minimizer.acc_factor = 1.001
	Jumol.run_cg(Minimizer, max_num_steps, eps)
	#
	# get the initial stress
	Jumol.calc_stress_tensor(molcalc)
	println("stress tensor:")
	pretty_table(molcalc.stress_tensor, noheader = true,
		         crop = :horizontal, formatters = ft_round(8))
	stress_before = molcalc.stress_tensor
	write(stress_strain_output_file, string(stress_before[1,1])*" "*string(stress_before[2,2])*" "*string(stress_before[1,2])*" ")
	#
	deform = Jumol.Aff_deform(ls.molstruc)
	lx0 = ls.molstruc.box.lx
	ly0 = ls.molstruc.box.ly
	delta_x = 0.1
	delta_stress = 1.0
	cntr_step = 1
	Borelax = Jumol.Box_relax(molcalc,deform,100,eps)
    Borelax.stress_tol = 1.0e-8
	while delta_stress > 0.0
		#
		# step number
		println("####")
		println("number of AQS step: ", cntr_step)
		# box deformation in transversal direction
		center_y = ls.molstruc.box.ly
		center_x = ls.molstruc.box.lx
		# deform box
		Jumol.set_tension_to_atom_group_2d(deform,center_x,center_y,delta_x,0.0)
		Jumol.relax_direction(Borelax, 2, 1.0e-2)
		# minimize potential energy
		Minimizer.alpha_min = 1.0e0
		Minimizer.acc_factor = 1.001
		Jumol.run_cg(Minimizer, max_num_steps, eps)
		# plot sample
		plot_sample(ls, ls.molstruc,
		            [lx0*(-0.05), lx0*(1.50), lx0*(1.50), lx0*(-0.05)],
		            [ly0*(-0.05), ly0*(-0.05), ly0*(1.05), ly0*(1.05)],
					7.0, 12.0)
		# get the stress
		Jumol.calc_stress_tensor(molcalc)
		println("stress tensor:")
	    pretty_table(molcalc.stress_tensor, noheader = true,
	        crop = :horizontal, formatters = ft_round(8))
		delta_stress = molcalc.stress_tensor[1,1] - stress_before[1,1]
		# update stress tensor
		stress_before = molcalc.stress_tensor
		#
		if delta_stress < 0
			write(stress_strain_output_file, string(cntr_step)*" "*string(stress_before[1,1])*" "*string(stress_before[2,2])*" "*string(stress_before[1,2])*"\n")
		end
		#println(delta_stress)
		#readline()
		cntr_step += 1
	end
	#if ls.save_def_vid
	#	gif(anim, video_filename*".gif", fps = 10)
	#end
end

## simple tension deformation with lammps
function run_simple_tension_deformation_lmp(ls::Local_structure, lmp_executable::String, stress_strain_output_file, loc_file_identification)
	# get the number of atoms in r_free
	cntr_atom_free = 0
	for atom in ls.molstruc.atom_list
		if atom.group == 0
			cntr_atom_free += 1
		end
	end
	#
	println("########")
	println("Local point scan:")
	println("num_x: ", ls.cntr_x, ", num_y: ", ls.cntr_y)
	println("coord_x: ", ls.x_coord, ", coord_y: ", ls.y_coord)
	println("theta: ", ls.theta)
	println("number of total atoms: ", length(ls.molstruc.atom_list))
	println("number of free atoms: ", cntr_atom_free)
	println("##")
	loc_sample_identification = "_cntrx"*string(ls.cntr_x)*"_cntry"*string(ls.cntr_y)*"_theta"*string(ls.theta)
	# write lammps input files
	Jumol.open_file(ls.molstruc.reader, "loc_sample_group_0.twodsilica", "w")
	Jumol.write_lmp_box_group(ls.molstruc.reader, ls.molstruc, 1, 0)
	Jumol.close_file(ls.molstruc.reader)
	Jumol.open_file(ls.molstruc.reader, "loc_sample_group_1.twodsilica", "w")
	Jumol.write_lmp_box_group(ls.molstruc.reader, ls.molstruc, 1, 1)
	Jumol.close_file(ls.molstruc.reader)
	# run the deformation with lammps
	run(`$lmp_executable -in deform.in`)
	println("lammps calculation finished")
	# read lammps stress output file
	strain_xx_vec, strain_yy_vec, stress_xx_vec, stress_yy_vec, stress_xy_vec, step_vec = read_lmp_stress_file(ls)
	# stress strain output file
	write(stress_strain_output_file, "#####\n")
	write(stress_strain_output_file, "number of atoms:\n")
	write(stress_strain_output_file, string(length(ls.molstruc.atom_list))*"\n")
	write(stress_strain_output_file, "number of atoms free:\n")
	write(stress_strain_output_file, string(cntr_atom_free)*"\n")
	write(stress_strain_output_file, "num_x_grid, num_y_grid, theta\n")
	write(stress_strain_output_file, string(ls.cntr_x)*" "*string(ls.cntr_y)*" "*string(ls.theta)*"\n")
	write(stress_strain_output_file, "pos_x, pos_y\n")
	write(stress_strain_output_file, string(ls.x_coord)*" "*string(ls.y_coord)*"\n")
	write(stress_strain_output_file, "# init_sigma_xx, init_sigma_yy, init_sigma_xy, sigma_xx, sigma_yy, sigma_xy\n")
	write(stress_strain_output_file, string(stress_xx_vec[1])*" "*string(stress_yy_vec[1])*" "*string(stress_xy_vec[1])*" "*
	                                 string(stress_xx_vec[end-1])*" "*string(stress_yy_vec[end-1])*" "*string(stress_xy_vec[end-1])*"\n")
	# video output filename
	video_output_filename = ls.output_path*"locdef_"*loc_file_identification*loc_sample_identification
	if ls.save_def_vid
		save_deformation_vid(ls, video_output_filename)
	end
end

## read lammps stress output
function read_lmp_stress_file(ls::Local_structure)
	stress_strain_output = Vector()
	stress_strain_file = open("stress_strain.dat","r")
	lines = readlines(stress_strain_file)
	close(stress_strain_file)
	#
	step_vec = Vector()
	strain_xx_vec = Vector(); strain_yy_vec = Vector()
	stress_xx_vec = Vector(); stress_yy_vec = Vector()
	stress_xy_vec = Vector()
	# read every second line only
	cntr = 2
	while cntr <= length(lines)
		data_vec = readdlm(IOBuffer(lines[cntr]))
		push!(strain_xx_vec, data_vec[1])
		push!(strain_yy_vec, data_vec[2])
		push!(stress_xx_vec, data_vec[3])
		push!(stress_yy_vec, data_vec[4])
		push!(stress_xy_vec, data_vec[5])
		push!(step_vec, Int(cntr/2))
		cntr += 2
	end
	return strain_xx_vec, strain_yy_vec, stress_xx_vec, stress_yy_vec, stress_xy_vec, step_vec
end

## save the deformation video
function save_deformation_vid(ls, filename)
	#load resopnse history
	video_molstruc = Jumol.Structure()
	video_molstruc.rc=2.5 # pseudo cutoff
	video_molstruc.rskin=0.1 # pseudo cutoff skin
	video_molstruc.pbx=0; video_molstruc.pby=0;
	## read lammpstrj
	initialize_structure_objects(video_molstruc)
	get_lines(video_molstruc.reader, "deform.lammpstrj")
	read_lammpstrj(video_molstruc.reader, video_molstruc, video_molstruc.reader.line_start_list[1])
	set_box_basis_vectors(video_molstruc.box)
	initialize_structure(video_molstruc, false)
	lx0 = video_molstruc.box.lx
	ly0 = video_molstruc.box.ly
	num_steps = length(video_molstruc.reader.line_start_list)
	anim = @animate for i in 1:num_steps
		read_lammpstrj(video_molstruc,"deform.lammpstrj",i)
		#mol_modifier = Modifier(video_molstruc)
		#define_group_out_of_circle(mol_modifier, video_molstruc.box.lx*0.5, video_molstruc.box.ly*0.5, ls.r_free)
		initialize_structure(video_molstruc, false)
		plot_sample(ls, video_molstruc,
		            [lx0*(-0.05), lx0*(1.50), lx0*(1.50), lx0*(-0.05)],
		            [ly0*(-0.05), ly0*(-0.05), ly0*(1.05), ly0*(1.05)])
	end
	gif(anim, filename*".gif", fps = 10)
end



## plot a sample
function plot_sample(ls::Local_structure, molstruc,
	                 border_x::Vector, border_y::Vector,
	                 size1=5.5, size2=3.5, lw=0.25,
	                 hfig=500, bfig=500,
					 col1="blue", col2="red")
    Plotter = Vis2d(molstruc)
    Plotter.hfig = hfig
    Plotter.bfig = bfig
    fig = plot_box(Plotter,"gray",1.0)
    plot_atomic_structure_binary(Plotter, size1, size2,
	                             border_x, border_y, col1, col2,
								 lw)
	# if frame_correction
	# 	plot!(Shape([molstruc.box.lx*(-0.04), molstruc.box.lx*(-0.0), molstruc.box.lx*(-0.0), molstruc.box.lx*(-0.04)],
	#     			[molstruc.box.ly*(-0.04), molstruc.box.ly*(-0.04), molstruc.box.ly*(1.04), molstruc.box.ly*(1.04)]),
	#     			color="white", linecolor="white", linewidth=0.0)
	#     plot!(Shape([molstruc.box.lx*(1.0), molstruc.box.lx*(1.04), molstruc.box.lx*(1.04), molstruc.box.lx*(1.0)],
	#     			[molstruc.box.ly*(-0.04), molstruc.box.ly*(-0.04), molstruc.box.ly*(1.04), molstruc.box.ly*(1.04)]),
	#     			color="white", linecolor="white", linewidth=0.0)
	#     plot!(Shape([molstruc.box.lx*(-0.04), molstruc.box.lx*(1.04), molstruc.box.lx*(1.04), molstruc.box.lx*(-0.04)],
	#     			[molstruc.box.ly*(-0.04), molstruc.box.ly*(-0.04), molstruc.box.ly*(-0.0), molstruc.box.ly*(-0.0)]),
	#     			color="white", linecolor="white", linewidth=0.0)
	#     plot!(Shape([molstruc.box.lx*(-0.04), molstruc.box.lx*(1.04), molstruc.box.lx*(1.04), molstruc.box.lx*(-0.04)],
	#     			[molstruc.box.ly*(1.0), molstruc.box.ly*(1.0), molstruc.box.ly*(1.04), molstruc.box.ly*(1.04)]),
	#     			color="white", linecolor="white", linewidth=0.0)
	#     Jumol.plot_box_border(Plotter, "black", 1.0)
	# end
    display(fig)
end
