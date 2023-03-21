## read LYS results
######################################
######################################
######################################
#
#
# (c) Franz Bamer, Ivy Wu Mai-2022
######################################


mutable struct LYS_Read
    filename::String
	num_x::Int64
	num_y::Int64
	num_theta::Int64
    LYS_Read(filename, num_x, num_y, num_theta) = new(filename, num_x, num_y, num_theta)
	noa_list::Vector
	noa_list_free::Vector
	num_matrix_angle_list::Vector
	pos_list::Vector
	init_stress_tensor_list::Vector #sigma_11, sigma_22, sigma_12
	fracture_stress_tensor_list::Vector #sigma_11, sigma_22, sigma_12
	delta_stress_tensor_list::Vector #sigma_11, simga_22, sigma_12
	#
	delta_stress_proj_list_normalized::Vector #sigma_11
	proj_delta_stress_list::Vector #sigma_11
	proj_min_number_list::Vector
	#
	red_factor_list::Vector # red factor for RGB plot
	blue_factor_list::Vector # blue factor for RGB plot
	#
	fig
	#
	sorted_soft_spots
	sorted_soft_spots_number_list
	CDF_soft_spots
	goodness
end

## read the output file
function read_lys_file(lr::LYS_Read)
	#
	println("####")
	println("reading file: "*lr.filename)
	println("num_x, num_y, num_theta")
	println(lr.num_x, ", ", lr.num_y, ", ", lr.num_theta)
	println("####")
	#
    file = open(lr.filename,"r")
	lines = readlines(file)
	close(file)
	#
	lr.noa_list = Vector()
	lr.noa_list_free = Vector()
	lr.num_matrix_angle_list = Vector()
	lr.pos_list = Vector()
	lr.init_stress_tensor_list = Vector()
	lr.fracture_stress_tensor_list = Vector()
	lr.delta_stress_tensor_list = Vector()
	cntr = 1
	for i in 1:lr.num_x
		for ii in 1:lr.num_y
			for iii in 1:lr.num_theta
				# empty line
				cntr += 1
				# number of atoms (description)
				cntr += 1
				# number of atoms (integer)
				push!(lr.noa_list, parse(Int, lines[cntr], base=10))
				cntr += 1
				# number of atoms free (description)
				cntr += 1
				# number of atoms free (integer)
				push!(lr.noa_list_free, parse(Int, lines[cntr], base=10))
				cntr += 1
				# num_x_grid, num_y_grid, theta (description)
				cntr += 1
				# num_x_grid, num_y_grid, theta (integer, integer, float)
				line_vec = readdlm(IOBuffer(lines[cntr]))
				nx = Int(line_vec[1])
				ny = Int(line_vec[2])
				angle = line_vec[3]
				push!(lr.num_matrix_angle_list, [nx, ny, angle])
				cntr += 1
				# pos_x, pos_y (description)
				cntr += 1
				# pos_x, pos_y (float, float)
				line_vec = readdlm(IOBuffer(lines[cntr]))
				pos_x = line_vec[1]
				pos_y = line_vec[2]
				push!(lr.pos_list,[pos_x, pos_y])
				cntr += 1
				# init_sigma_xx, init_simga_yy, init_simga_xy, sigma_xx, sigma_yy, sigma_xy (description)
				cntr += 1
				# init_sigma_xx, init_simga_yy, init_simga_xy, sigma_xx, sigma_yy, sigma_xy (floats)
				line_vec = readdlm(IOBuffer(lines[cntr]))
				init_sigma_xx = line_vec[1]
				init_sigma_yy = line_vec[2]
				init_sigma_xy = line_vec[3]
				push!(lr.init_stress_tensor_list, [init_sigma_xx, init_sigma_yy, init_sigma_xy])
				sigma_xx = line_vec[4]
				sigma_yy = line_vec[5]
				sigma_xy = line_vec[6]
				push!(lr.fracture_stress_tensor_list, [sigma_xx, sigma_yy, sigma_xy])
				cntr += 1
				#
				push!(lr.delta_stress_tensor_list, lr.fracture_stress_tensor_list[end] - lr.init_stress_tensor_list[end])
			end
		end
	end
end

## project yield stress into the main principle shear stress directions
function project_delta_stress(lr::LYS_Read)
	lr.proj_delta_stress_list = Vector()
	lr.proj_min_number_list = Vector()
	cntr = 1
	theta_l = 0.0 #FIXME correct this later so that it can be changed
	for i in 1:lr.num_x
		for ii in 1:lr.num_y
			##
			# first value -> create variable of proj_delta_stress
			min_proj_delta_stress = 1.0e10
			proj_number = 0
			for iii in 1:lr.num_theta
				sig_11 = lr.delta_stress_tensor_list[cntr][1]
				sig_22 = lr.delta_stress_tensor_list[cntr][2]
				sig_12 = lr.delta_stress_tensor_list[cntr][3]
				shear_stress = (sig_11 - sig_22)*0.5 #FIXME
				theta = lr.num_matrix_angle_list[cntr][3]
				alpha = theta_l - theta
				if abs(alpha) == 0.0
					proj_delta_stress = sig_11 / cosd(alpha)
					min_proj_delta_stress = proj_delta_stress
					proj_number = cntr
					# remaining values for theta
					if proj_delta_stress < min_proj_delta_stress
						min_proj_delta_stress = proj_delta_stress
						proj_number = cntr
					end
				end
				cntr += 1
			end
			push!(lr.proj_delta_stress_list, min_proj_delta_stress)
			push!(lr.proj_min_number_list, proj_number)
		end
	end
end

## get colorscheme
function plot_lys_res(lr::LYS_Read, r_free::Float64, num_sample::Int64, sig::Int64, molstruc)
	## get maximum delta tensile stress
	max_delta_sig_proj = -1.0e10
	for delta_sig_proj in lr.proj_delta_stress_list
		if delta_sig_proj > max_delta_sig_proj
			max_delta_sig_proj = delta_sig_proj
		end
	end
	## get minimum delta tensilr stress
	min_delta_sig_proj = 1.0e10
	for delta_sig_proj in lr.proj_delta_stress_list
		if delta_sig_proj < min_delta_sig_proj
			min_delta_sig_proj = delta_sig_proj
		end
	end
	#
	lr.delta_stress_proj_list_normalized = Vector()
	lr.red_factor_list = Vector()
	lr.blue_factor_list = Vector()
	for dsig in lr.proj_delta_stress_list
		# normalize between 0 and 1
		norm_val = (dsig[1] - min_delta_sig_proj) / (max_delta_sig_proj - min_delta_sig_proj)
		push!(lr.delta_stress_proj_list_normalized, norm_val)
		push!(lr.red_factor_list, 1 - norm_val)
		push!(lr.blue_factor_list, norm_val)
	end
	# plot the points in the sample
	lr.fig = plot(border=:none,aspect_ratio=1,
	            legend=false,fmt=:pdf, title="r_free = "*string(r_free)*", num_sample = "*string(num_sample)*", sig = "*string(sig))
	plot!(size=(1000,1000))
	num_points = lr.num_x*lr.num_y
	for i in 1:num_points
		num = lr.proj_min_number_list[i]
		x = lr.pos_list[num][1]
		y = lr.pos_list[num][2]
		fr = lr.red_factor_list[i]
		fb = lr.blue_factor_list[i]
		scatter!([x],[y], markersize=13.0, markershape=:square,
		         color=RGB(fr, 0.0, fb), markerstrokecolor=RGB(fr, 0.0, fb))
		## plot periodic result if applicable
 		if lr.num_matrix_angle_list[num][1] == 1
 			xp = x + molstruc.box.lx
			scatter!([xp],[y], markersize=13.0, markershape=:square,
			         color=RGB(fr, 0.0, fb), markerstrokecolor=RGB(fr, 0.0, fb))
		end
 		if lr.num_matrix_angle_list[num][2] == 1
 			yp = y + molstruc.box.ly
			scatter!([x],[yp], markersize=13.0, markershape=:square,
			         color=RGB(fr, 0.0, fb), markerstrokecolor=RGB(fr, 0.0, fb))
 		end
		if lr.num_matrix_angle_list[num][1] == 1 && lr.num_matrix_angle_list[num][2] == 1
			xp = x + molstruc.box.lx
			yp = y + molstruc.box.ly
			scatter!([xp],[yp], markersize=13.0, markershape=:square,
			         color=RGB(fr, 0.0, fb), markerstrokecolor=RGB(fr, 0.0, fb))
		end
	end
end

## find the first soft spots
function sort_soft_spots(lr::LYS_Read)
	#
	permut = sortperm(lr.delta_stress_proj_list_normalized)
	# sort the soft spot lists
	lr.sorted_soft_spots = lr.delta_stress_proj_list_normalized[permut]
	lr.sorted_soft_spots_number_list = lr.proj_min_number_list[permut]
end

## calculate prediction goodness
function calc_prediction_goodness(lr::LYS_Read, event_coords, num_events)
	if num_events > size(event_coords,1)
		num_events = size(event_coords,1)
	end
	lr.goodness = zeros(num_events,2)
	for i in 1:num_events
		#
		lr.goodness[i,1] = i
		# event position
		event_pos = event_coords[i,:]
		# find nearest grid point
		num_points = lr.num_x*lr.num_y
		min_dist = 99999999.99
		# start searching values
		event_grid_point = lr.pos_list[1]
		local_event_stress = lr.delta_stress_proj_list_normalized[1]
		for ii in 1:num_points
			num = lr.proj_min_number_list[ii]
			grid_pos = lr.pos_list[num]
			dist = norm(event_pos - grid_pos)
			if dist < min_dist
				event_grid_point = copy(grid_pos)
				local_event_stress = lr.delta_stress_proj_list_normalized[ii]
				min_dist = dist
			end
		end
		#
		#println(event_grid_point)
		#println(local_event_stress)
		#readline()
		lr.goodness[i,2] = 1.0 - 2*call_CDF(lr, local_event_stress)
	end

end

## plot the soft spots
function plot_soft_spots(lr::LYS_Read, num, dx::Float64=5.0, dy::Float64=-1.0)
	for i in 1:num
		coord_index = lr.sorted_soft_spots_number_list[i]
		x = lr.pos_list[coord_index][1]
		y = lr.pos_list[coord_index][2]
		scatter!([x], [y], markershape=:star4, markersize=12.0,
		                   color="yellow", markerstrokecolor="yellow", markerstrokewidth=2.0)
		#annotate!([x+dx], [y+dy], text(string(i), color="green"))
	end
end

## load and plot the aposteriori extracted event spots
function plot_event_spots(lr::LYS_Read, event_coords, num, dx::Float64=5.0, dy::Float64=-1.0, size=5.0)
	for i in 1:num
		x = event_coords[i,1]
		y = event_coords[i,2]
		scatter!([x], [y], markershape=:circle, markersize = size,
		                   color="green", markerstrokecolor="white", markerstrokewidth=0.5)
		annotate!([x+dx], [y+dy], text(string(i), color="black"))
	end
end

## plot the molecular structure
function plot_molecular_structure(lr::LYS_Read, molstruc, PI, op=1.0)
	Plotter = Vis2d(molstruc)
	Plotter.hfig=PI.hfig; Plotter.bfig=PI.bfig
	lx = molstruc.box.lx
	ly = molstruc.box.ly
	plot_atomic_structure_binary(Plotter, PI.at_size1, PI.at_size2,
	                             PI.factor_vec_x*lx,
								 PI.factor_vec_y*ly,
								 "lightblue", "lightgray",
								 0.25, op)
	## cut of unnecessary frame
	xl1 = 0.0 - lx*0.03
	x1 = 0.0
	x2 = lx
	xr2 = lx*1.03
	yb1 = 0.0 - ly*0.03
	y1 = 0.0
	y2 = ly
	yt2 = ly*1.03
	plot!(Shape([xl1,x1,x1,xl1],[yb1,yb1,yt2,yt2]), color="white", linewidth=0.0, linecolor="white")
	plot!(Shape([xl1,xr2,xr2,xl1],[yb1,yb1,y1,y1]), color="white", linewidth=0.0, linecolor="white")
	plot!(Shape([x2,xr2,xr2,x2],[yb1,yb1,yt2,yt2]), color="white", linewidth=0.0, linecolor="white")
	plot!(Shape([xl1,xr2,xr2,xl1],[y2,y2,yt2,yt2]), color="white", linewidth=0.0, linecolor="white")
	## plot box border
	plot!(Shape([x1,x2,x2,x1],[y1,y1,y2,y2]), linewidth=3.0, opacity=0.0)
end

## display the figure
function display_fig(lr::LYS_Read, sig, num_sample, r_free, save_figure::Bool=false)
	if save_figure
		savefig(lr.fig, "vis_output/lys_res_sig"*string(sig)*"_num"*string(num_sample)*"_r"*string(r_free)*".pdf")
	end
	display(lr.fig)
end

## calculate the CDF of soft spots
function calc_CDF(lr::LYS_Read)
	# sort local yield stresses
	num = lr.num_x*lr.num_y
	lr.CDF_soft_spots = zeros(num,2)
	for i in 1:num
		lr.CDF_soft_spots[i,1] = lr.sorted_soft_spots[i]
		lr.CDF_soft_spots[i,2] = float(i)/float(num)
	end
end

## look-up value from the CDF
function call_CDF(lr::LYS_Read, value)
	if value < lr.CDF_soft_spots[1,1]
		return 0.0
	end
	num = lr.num_x*lr.num_y
	for i in 1:num-1
		low_val = lr.CDF_soft_spots[i,1]
		high_val = lr.CDF_soft_spots[i+1,1]
		if value >= low_val && value <= high_val
			return (lr.CDF_soft_spots[i,2]+lr.CDF_soft_spots[i+1,2])*0.5
		end
	end
	if value > lr.CDF_soft_spots[end,1]
		return 1.0
	end
end
