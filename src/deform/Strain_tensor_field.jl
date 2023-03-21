## Strain tensor field class
######################################
######################################
######################################
#
#
# (c) Franz Bamer Sep-2022
######################################

using LinearAlgebra
#using VegaLite, DataFrames

mutable struct Strain_tensor_field
    molstruc_step1
    molstruc_step2
    Strain_tensor_field(molstruc_step1,
               molstruc_step2) = new(molstruc_step1,
                                      molstruc_step2)
	loc_def_grad_list::Vector
	loc_Right_Cauchy_Green_list::Vector
	loc_Green_Lag_strain_list::Vector
	loc_Green_Lag_principal_strain_list::Vector
	loc_Green_Lag_principal_direc_list::Vector
	loc_Green_Lag_extr_shear_strain_list::Vector
	#
	loc_Left_Cauchy_green_list::Vector
	loc_Alm_strain_list::Vector
end

## get affine displacement field from the boxes
# deformation gradient is saved as a list of 2x2 matrices
function calc_def_grad_field(sf::Strain_tensor_field)
	sf.loc_def_grad_list = Vector()
	for i in 1:sf.molstruc_step1.noa
		atom_center_ref = sf.molstruc_step1.atom_list[i]
		center_num = atom_center_ref.number
		#
		dist_mat_ref = zeros(length(atom_center_ref.neighbor_indices),2)
		dist_mat_cur = zeros(length(atom_center_ref.neighbor_indices),2)
		for ii in 1:length(atom_center_ref.neighbor_indices)
			neigh_num = atom_center_ref.neighbor_indices[ii]
			dist_vec_ref = get_distance_between2atoms(sf.molstruc_step1, center_num, neigh_num)
			dist_vec_cur = get_distance_between2atoms(sf.molstruc_step2, center_num, neigh_num)
			#
			dist_mat_ref[ii,:] = dist_vec_ref[1:2]
			dist_mat_cur[ii,:] = dist_vec_cur[1:2]
		end
		push!(sf.loc_def_grad_list,inv(transpose(dist_mat_ref)*dist_mat_ref)*transpose(dist_mat_ref)*dist_mat_cur)
	end
end

## get affine displacement field from the boxes
# deformation gradient is saved as a list of 2x2 matrices
function calc_Green_Lagrange_field(sf::Strain_tensor_field)
	sf.loc_Right_Cauchy_Green_list = Vector()
	sf.loc_Green_Lag_strain_list = Vector()
	sf.loc_Green_Lag_principal_strain_list = Vector()
	sf.loc_Green_Lag_principal_direc_list = Vector()
	sf.loc_Green_Lag_extr_shear_strain_list = Vector()
	I2 = Matrix{Float64}(I, 2, 2)
	for loc_F in sf.loc_def_grad_list
		# Right Cauchy Green
		C = transpose(loc_F)*loc_F
		push!(sf.loc_Right_Cauchy_Green_list, C)
		# Green Lagrange
		E = 0.5*(C-I2)
		push!(sf.loc_Green_Lag_strain_list, E)
		# principal strains (E is symmetric)
		first_val = (E[1,1]-E[2,2])*0.5
		sqrtval = sqrt( (E[1,1]-E[2,2])*(E[1,1]-E[2,2])*0.25 + E[1,2]*E[1,2] )
		eps1 = first_val + sqrtval
		eps2 = first_val - sqrtval
		eps_vec = zeros(2)
		eps_vec[1] = eps1
		eps_vec[2] = eps2
		push!(sf.loc_Green_Lag_principal_strain_list, eps_vec)
		# principal strain axes (E is symmetric)
		leng1 = sqrt(1+E[1,2]*E[1,2]/(E[2,2]-eps1)^2)
		lambda_mat = zeros(2,2)
		lambda_mat[1,1] = 1.0/leng1
		lambda_mat[2,1] = -E[1,2]/(E[2,2]-eps1)/leng1
		leng2 = sqrt(1+E[1,2]*E[1,2]/(E[2,2]-eps2)^2)
		lambda_mat[1,2] = 1.0/leng2
		lambda_mat[2,2] = -E[1,2]/(E[2,2]-eps2)/leng2
		push!(sf.loc_Green_Lag_principal_direc_list,lambda_mat)
		#
		extr_shear_strain_mat = zeros(1,1)
		extr_shear_strain_mat[1,1] = (eps2 - eps1)*0.5
		push!(sf.loc_Green_Lag_extr_shear_strain_list, extr_shear_strain_mat)
	end
end


## plot the unit box
function plot_unit_box(sf::Strain_tensor_field, linewidth=0.5)
    p_mat = zeros(3,5)
	p_mat[:,2] = zeros(3) + [1;0;0]
	p_mat[:,3] = p_mat[:,2] + [0;1;0]
	p_mat[:,4] = p_mat[:,3] - [1;0;0]
	unit_box_plot = plot(p_mat[1,:],p_mat[2,:],border=:none,aspect_ratio=1,
	            legend=false,color="black",lw=linewidth,fmt=:pdf)
	plot!(size=(1000,1000))
	return unit_box_plot
end

## plot a displacement field
# create a plotting object before calling the function
function plot_field(sf::Strain_tensor_field, field::Vector, entries::Vector{Int64},
	                boxlw::Float64)
	# normalize strain vector
	strain_vec = zeros(length(field))
	for i in 1:length(field)
		strain_vec[i] = field[i][entries[1],entries[2]]
	end
	norm_strain_vec, red_vec, blue_vec, norm_int = normalize_strain_vec(sf, strain_vec)
	# normalize box coordinates
	H_1 = zeros(3,3)
    H_1[:,1] = sf.molstruc_step1.box.h1
    H_1[:,2] = sf.molstruc_step1.box.h2
    H_1[:,3] = sf.molstruc_step1.box.h3
	# map the coordinates from step 1 to the unit box
	coord_mat_unit_square = zeros(sf.molstruc_step1.noa,3)
	H_1_inv = inv(H_1)
    for i in 1:sf.molstruc_step1.noa
        coord_mat_unit_square[i,:] = H_1_inv*sf.molstruc_step1.atom_list[i].pos
    end
	# plot the field
	cm = colormap("RdBu",1001,logscale=false)
	fig = plot_unit_box(sf, boxlw)
	size1 = 11.0
	size2 = 21.0
	for i in 1:size(coord_mat_unit_square,1)
		atom_type = sf.molstruc_step1.atom_list[i].type
		size = 0
		if atom_type  == 1
			size = size1
		end
		if atom_type == 2
			size = size2
		end
		x = coord_mat_unit_square[i,1]
		y = coord_mat_unit_square[i,2]
		fr = red_vec[i]
		fb = blue_vec[i]
		intcolval = norm_int[i]
		scatter!([x],[y], markersize=size, markershape=:circle,
		         color=cm[intcolval], markerstrokecolor=cm[intcolval])
	end
	return fig
end





## normalize strain values
function normalize_strain_vec(sf::Strain_tensor_field, strain_vec)
	## get minimum strain value
	min_strain = 1.0e10
	for strain in strain_vec
		if strain < min_strain
			min_strain = strain
		end
	end
	## get maximum strain value
	max_strain = -1.0e10
	for strain in strain_vec
		if strain > max_strain
			max_strain = strain
		end
	end
	#
	strain_list_normalized = Vector()
	red_factor_list = Vector()
	blue_factor_list = Vector()
	strain_int = Vector()
	for strain in strain_vec
		# normalize between 0 and 1
		norm_val = (strain - min_strain) / (max_strain - min_strain)
		push!(strain_list_normalized, norm_val)
		push!(red_factor_list, 1 - norm_val)
		push!(blue_factor_list, norm_val)
		push!(strain_int, floor(Int,norm_val*1000.0+1))
	end
	return strain_list_normalized, red_factor_list, blue_factor_list, strain_int
end


## normalize strain values around zero
function normalize_strain_vec_around_zero(sf::Strain_tensor_field, strain_vec)
	#
	strain_list_normalized = Vector()
	red_factor_list = Vector()
	blue_factor_list = Vector()
	strain_int = Vector()
	## get minimum strain value
	min_strain = 1.0e10
	for strain in strain_vec
		if strain < min_strain
			min_strain = strain
		end
	end
	for strain in strain_vec
		if strain < 0
			# normalize between 0 and 0.5 for all negative strain values
			norm_val = (strain - min_strain) / (0 - min_strain)
			push!(strain_list_normalized, norm_val)
			push!(red_factor_list, 1 - norm_val)
			push!(blue_factor_list, norm_val)
			push!(strain_int, floor(Int,norm_val*1000.0+1))
		end
	end
	## get maximum strain value
	max_strain = -1.0e10
	for strain in strain_vec
		if strain > max_strain
			max_strain = strain
		end
	end
	for strain in strain_vec
		if strain > 0
			# normalize between 0 and 0.5 for all negative strain values
			norm_val = (strain - 0) / (max_strain - 0)
			push!(strain_list_normalized, norm_val)
			push!(red_factor_list, 1 - norm_val)
			push!(blue_factor_list, norm_val)
			push!(strain_int, floor(Int,norm_val*1000.0+1))
		end
	end
	return strain_list_normalized, red_factor_list, blue_factor_list, strain_int
end
