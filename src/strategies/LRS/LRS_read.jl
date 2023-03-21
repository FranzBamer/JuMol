## Read LRS results from file
######################################
######################################
######################################
# read results
# investigate directional data
#
# (c) Franz Bamer, Ivy Wu September-2022
######################################


using LinearAlgebra

## LRS reader
mutable struct LRS_Read
	fname::String
	num_x::Int64
	num_y::Int64
    LRS_Read(fname, num_x, num_y) = new(fname, num_x, num_y)
	grid_num_list::Vector
	pos_list::Vector
	mu_list::Vector
	sig_list::Vector
	physical_nodes_mat_list::Vector
	bond_mat_list::Vector
	bond_list_list::Vector
	bond_n_list_list::Vector
	bond_l_list_list::Vector
	hist_list::Vector
	num_intersec_hist::Int64
	pdf_hist::Vector
	ring_Si_list_list_list::Vector
	#
	CDF_het::Array
	#
	F_mean_list::Vector
	F_mean_n_list::Vector
	F_sig_list::Vector
	F_sig_n_list::Vector
	F_mean_diff_list::Vector
	F_sig_diff_list::Vector
	lambda_mean_list::Vector
	eigmat_mean_list::Vector
	lambda_sig_list::Vector
	eigmat_sig_list::Vector
	#
	F_mean_ring_list_list::Vector
	F_sig_ring_list_list::Vector
end

function read_file(lr::LRS_Read)
	file = open(lr.fname, "r")
	lines = readlines(file)
	close(file)
	#
	lr.grid_num_list = Vector()
	lr.pos_list = Vector()
	lr.mu_list = Vector()
	lr.sig_list = Vector()
	lr.physical_nodes_mat_list = Vector()
	lr.bond_mat_list = Vector()
	lr.ring_Si_list_list_list = Vector()
	#
	cntr = 1
	for i in 1:lr.num_x
		for ii in 1:lr.num_y
			# num_x_grid, num_y_grid
			cntr += 1
			# read
			line_vec = readdlm(IOBuffer(lines[cntr]))
			nx = Int(line_vec[1])
			ny = Int(line_vec[2])
			push!(lr.grid_num_list, [nx, ny])
			cntr += 1
			#println("num_x_grid, num_y_grid")
			#println(string(nx)*", "*string(ny))
			# pos_x, pos_y
			cntr += 1
			# read
			line_vec = readdlm(IOBuffer(lines[cntr]))
			pos_x = line_vec[1]
			pos_y = line_vec[2]
			push!(lr.pos_list,[pos_x, pos_y])
			cntr += 1
			#println("pos_x, pos_y")
			#println(string(pos_x)*", "*string(pos_y))
			# ring statistics mu, sigma
			cntr += 1
			# read
			line_vec = readdlm(IOBuffer(lines[cntr]))
			push!(lr.mu_list, line_vec[1])
			push!(lr.sig_list, line_vec[2])
			cntr += 1
			#println("mu, sig")
			#println(line_vec)
			# physical nodes in matrix form: line1->x_coords, line2->y-coords
			cntr += 1
			# read
			physical_nodes_mat_x = readdlm(IOBuffer(lines[cntr]))
			physical_nodes_mat = zeros(2, length(physical_nodes_mat_x))
			physical_nodes_mat[1,:] = physical_nodes_mat_x
			cntr += 1
			physical_nodes_mat_y = readdlm(IOBuffer(lines[cntr]))
			physical_nodes_mat[2,:] = physical_nodes_mat_y
			push!(lr.physical_nodes_mat_list, physical_nodes_mat)
			cntr += 1
			#println("physical nodes matrix")
			#println(physical_nodes_mat)
			# bond list in matrix form: line1->x_components, line2->y-components
			cntr += 1
			# read
			bond_mat_x = readdlm(IOBuffer(lines[cntr]))
			bond_mat = zeros(2, length(bond_mat_x))
			bond_mat[1,:] = bond_mat_x
			cntr += 1
			bond_mat_y = readdlm(IOBuffer(lines[cntr]))
			bond_mat[2,:] = bond_mat_y
			push!(lr.bond_mat_list, bond_mat)
			cntr += 1
			#println("bond list matrix")
			#println(bond_mat)
			# Si ring nodes in matrix form: line1->x_coords, line2->y-coords (2 lines for one ring)
			cntr += 1
			# read number
			num_rings = parse(Int, lines[cntr], base=10)
			#println("number of rings:, ", num_rings)
			cntr += 1
			# read Si atoms of the rings
			ring_Si_list_list = Vector()
			for j in 1:num_rings
				Si_x_coords = readdlm(IOBuffer(lines[cntr]))
				cntr += 1
				Si_y_coords = readdlm(IOBuffer(lines[cntr]))
				cntr += 1
				#println("Ring number: ", j)
				#println(Si_x_coords)
				#println(Si_y_coords)
				Si_pos_list = Vector()
				for jj in 1:length(Si_x_coords)
					x = Si_x_coords[jj]
					y = Si_y_coords[jj]
					push!(Si_pos_list,[x,y])
				end
				push!(ring_Si_list_list, Si_pos_list)
			end
			push!(lr.ring_Si_list_list_list, ring_Si_list_list)

			#readline()
		end
	end
end

## include plotting functions
include("LRS_plottings.jl")

## include general functions (mat to list, ect.)
include("LRS_general_functions.jl")

## include fabric tensor calculations
include("LRS_empirical_distribution_density.jl")


## get the list of bonds for every area
function extract_bond_lists(lr::LRS_Read)
	lr.bond_list_list = Vector()
	lr.bond_n_list_list = Vector()
	lr.bond_l_list_list = Vector()
	for mat in lr.bond_mat_list
		bond_list = mat2vec_list(mat)
		bond_n_list, bond_l_list = normalize_list(bond_list)
		for i in 1:length(bond_l_list)
			bond_l_list[i] = bond_l_list[i] / 3.05
		end
		# get the bonds in both directions for the symmetric distribution
		lbl = length(bond_list)
		for i in 1:lbl
			push!(bond_list, bond_list[i]*(-1.0))
			push!(bond_n_list, bond_n_list[i]*(-1.0))
			push!(bond_l_list, bond_l_list[i])
		end
		push!(lr.bond_list_list, bond_list)
		push!(lr.bond_n_list_list, bond_n_list)
		push!(lr.bond_l_list_list, bond_l_list)
	end
end

## extract the empirical distribution density
function extract_emp_dist_dens_list(lr::LRS_Read, num_intersec)
	lr.num_intersec_hist = num_intersec
	lr.hist_list = Vector()
	for bond_list in lr.bond_n_list_list
		push!(lr.hist_list, calc_empirical_distribution_density(lr.num_intersec_hist, bond_list))
	end
end


## predict anisotropy
function predict_anisotropic_bond_direc(lr::LRS_Read, type::String)
	#
	extract_bond_lists(lr)
	#
	lr.F_mean_list = calc_fabric_tensor_list(lr.bond_list_list)
	lr.F_mean_n_list = calc_fabric_tensor_list(lr.bond_n_list_list)
	lr.F_mean_diff_list = Vector()
	for i in 1:length(lr.F_mean_list)
		push!(lr.F_mean_diff_list, lr.F_mean_list[i] - lr.F_mean_n_list[i])
	end
	#
	lr.F_sig_list = calc_fabric_tensor_list_sig(lr.F_mean_list, lr.bond_list_list)
	lr.F_sig_n_list = calc_fabric_tensor_list_sig(lr.F_mean_n_list, lr.bond_list_list)
	lr.F_sig_diff_list = Vector()
	for i in 1:length(lr.F_sig_list)
		push!(lr.F_sig_diff_list, lr.F_sig_list[i] - lr.F_sig_n_list[i])
	end
	#
	alpha_list = Vector()
	#
	lr.lambda_mean_list = Vector()
	lr.eigmat_mean_list = Vector()
	lr.lambda_sig_list = Vector()
	lr.eigmat_sig_list = Vector()
	cntr = 1
	for i in 1:length(lr.F_mean_list)
		F_mean = lr.F_mean_list[cntr]
		F_sig = lr.F_sig_list[cntr]
		#
		lambdavec_mean, eigmat_mean = calc_eigensolution_two_times_two(F_mean)
		push!(lr.lambda_mean_list, lambdavec_mean)
		push!(lr.eigmat_mean_list, eigmat_mean)
		vec1_mean = eigmat_mean[:,1]
		vec2_mean = eigmat_mean[:,2]
		# transform standard deviation into the principal axes space
		lambdavec_sig, eigmat_sig = calc_eigensolution_two_times_two(F_sig)
		push!(lr.lambda_sig_list, lambdavec_sig)
		push!(lr.eigmat_sig_list, eigmat_sig)
		vec1_sig = eigmat_sig[:,1]
		vec2_sig = eigmat_sig[:,2]
		#
		e1 = zeros(2); e1[1] = 1.0; e1[2] = 0.0
		#
		#tanalpha = 1.0/(lambdavec[1] * cos(vec1[1]/norm(vec1)) ) * abs(atan(vec1[2] / vec1[1]))
		#alpha_mean = 1.0/(lambdavec_mean[1]) * abs(atan(vec1_mean[2] / vec1_mean[1]))
		if type == "mean"
			alpha = 1.0/abs(dot(vec1_mean, e1))
		elseif type == "sig"
			#alpha = abs(atan(vec1_sig[2] / vec1_sig[1]))
			alpha = 1.0/abs(dot(vec1_sig*lambdavec_sig[1], e1))
		elseif type == "mean_and_sig"
			alpha = 1.0/( abs(dot(vec1_sig*lambdavec_sig[1], vec1_mean*lambdavec_mean[1])) ) * 1.0/abs(dot(vec1_mean*lambdavec_mean[1], e1))
			#println( (1.0/abs(dot(vec1_mean, e1))) / (1.0/abs(dot(vec1_sig, e1))) )
			#readline()
		else
			println("not allowed")
		end
		push!(alpha_list, alpha)
		#
		cntr += 1
	end
	return alpha_list
end


## calculate the CDF of a vector
function calc_CDF(vector)
	permut = sortperm(vector)
	# sort the soft spot lists
	vector_sorted = vector[permut]
	# sort local yield stresses
	num = length(vector)
 	CDF_mat= zeros(num,2)
	for i in 1:num
		CDF_mat[i,1] = vector_sorted[i]
		CDF_mat[i,2] = float(i)/float(num)
	end
	return CDF_mat
end


## look-up value from the CDF
function call_CDF(value, CDF_mat)
	if value < CDF_mat[1,1]
		return 0.0
	end
	num = size(CDF_mat,1)
	for i in 1:num-1
		low_val = CDF_mat[i,1]
		high_val = CDF_mat[i+1,1]
		if value >= low_val && value <= high_val
			return (CDF_mat[i,2]+CDF_mat[i+1,2])*0.5
		end
	end
	if value > CDF_mat[end,1]
		return 1.0
	end
end


## calculate prediction goodness
function calc_prediction_goodness_event_to_grid(lr::LRS_Read, event_coords, num_events, pred_vector, CDF_mat)
	if num_events > size(event_coords,1)
		num_events = size(event_coords,1)
	end
	goodness = zeros(num_events,2)
	for i in 1:num_events
		#
		goodness[i,1] = i
		# event position
		event_pos = event_coords[i,:]
		# find nearest grid point
		num_points = lr.num_x*lr.num_y
		min_dist = 99999999.99
		# start searching values
		event_grid_point = lr.pos_list[1]
		local_val = pred_vector[1]
		for ii in 1:num_points
			grid_pos = lr.pos_list[ii]
			dist = norm(event_pos - grid_pos)
			if dist < min_dist
				event_grid_point = copy(grid_pos)
				local_val = pred_vector[ii]
				min_dist = dist
			end
		end
		#
		#println(event_grid_point)
		#println(local_event_stress)
		#readline()
		#println("local_val")
		#println(local_val)
		#println("CDF_mat")
		#println(CDF_mat)
		#readline()
		goodness[i,2] = 1.0 - 2*call_CDF(local_val, CDF_mat)
	end
	return goodness
end





####################################
####################################
#### ring-wise calculations
## predict mean anisotropy of every ring
function predict_mean_anisotropy_rings(lr::LRS_Read)
	#
	calc_fabric_ring_tensor_list(lr)
	#
	tan_alpha_list = Vector()
	cntr_region = 1
	for F_list in lr.F_mean_ring_list_list
		mean_tan_alpha = 0.0
		cntr_ring = 1
		for F in F_list
			size_ring = length(lr.ring_Si_list_list_list[cntr_region][cntr_ring])
			lambdavec, eigmat = calc_eigensolution_two_times_two(F)
			vec1 = eigmat[:,1]
			vec2 = eigmat[:,2]
			mean_tan_alpha += abs(atan(vec1[2] / vec1[1])) * 1.0/lambdavec[1] * 1.0/float(size_ring)
			#
			cntr_ring += 1
		end
		mean_tan_alpha = mean_tan_alpha / float(length(F_list))
		push!(tan_alpha_list, mean_tan_alpha)
		#
		cntr_region += 1
	end
	return tan_alpha_list
end


## calculate a list of fabric tensors for every ring
function calc_fabric_ring_tensor_list(lr::LRS_Read)
	lr.F_mean_ring_list_list = Vector()
	lr.F_sig_ring_list_list = Vector()
	for loc_region in lr.ring_Si_list_list_list
		F_mean_ring_list = Vector()
		F_sig_ring_list = Vector()
		for ring_list in loc_region
			bond_list = Vector()
			for i in 2:length(ring_list)
				push!(bond_list, ring_list[i] - ring_list[i-1])
			end
			push!(bond_list, ring_list[1] - ring_list[end])
			# evaluate fabric tensor of one ring
			F_mean_ring, F_sig_ring = calc_fabric_tensor(bond_list)
			push!(F_mean_ring_list, F_mean_ring)
			push!(F_sig_ring_list, F_sig_ring)
		end
		push!(lr.F_mean_ring_list_list, F_mean_ring_list)
		push!(lr.F_sig_ring_list_list, F_sig_ring_list)
	end
end
