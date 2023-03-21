## struct for LRS method (local ring statistics)
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
# (c) Franz Bamer, Ivy Wu September-2022
######################################

# this file must be generalized to local probing in general
# local potential energy, local free volume, local ring statistics,
# local ring neighborhood statistics, local fabric tensors,

include("../dual_structure/Dual_network.jl")
## lrs method
mutable struct LRS
	LI #(struct defined by input)
	CI #(struct defined by input)
	SI #(struct defined by input)
	PI #(struct defined by input)
	r_free::Float64
	num_sample::Int64
	target_params
	allowed_r_list
	target_stat_list
    LRS(LI, CI, SI, PI, r_free, num_sample,
	    target_params, allowed_r_list, target_stat_list) = new(LI, CI, SI, PI, r_free, num_sample,
															   target_params, allowed_r_list, target_stat_list)
	## global network
	nw
	## scanning parameters
	x_coord_list::Vector
	y_coord_list::Vector
end

## load the sample
function load_dualnetwork(lrs::LRS, visualize=false)
	file_load = lrs.SI.input_path*lrs.SI.filename_dualnet
	println("... importing sample ...")

	# size of structure
	bond_length= 3.05
	minimize_factor = 1
	# Aboav-Weaire parameter
	alpha4 = 0.3
	alpha5 = 0.3

	target_sig = lrs.target_params[lrs.SI.sig,2]
	r_list = lrs.allowed_r_list[lrs.SI.sig]
	target_stat = lrs.target_stat_list[lrs.SI.sig,:]
	alpha_list = ones(length(r_list))*alpha4

	nw = Network(bond_length,"silica")
	nw.full_min_rate_factor = minimize_factor
	nw.ring_size_list_allowed = r_list
	nw.ring_stat_vec = target_stat

	# initialize dual network
	Dualnetwork = Dual_network(nw.bond_length)
	read_dual_network(Dualnetwork, file_load)

	# create the physical equivalent
	Physnetwork = Phys_network(Dualnetwork)
	Physnetwork.include_lammps = lrs.CI.include_lmp
	Physnetwork.lammps_executable = lrs.CI.lmp_executable
	Physnetwork.sig = lrs.SI.sig

	Jumol.dual2phys(Physnetwork)
	Jumol.build_phys_molstruc(Physnetwork)
	Jumol.minimize_phys_network(Physnetwork)
	Jumol.update_node_pos(Physnetwork)
	Jumol.build_silica_molstruc(Physnetwork)
	coordination_check = Jumol.minimize_silica_network_lammps_plus_coord_check(Physnetwork)
	# update the final silica positions in all networks
	Jumol.update_silica_node_pos(Physnetwork)
	Jumol.build_phys_molstruc(Physnetwork)
	Jumol.phys2dual(Dualnetwork, Physnetwork)
	#
	nw.Dualnetwork = Dualnetwork
	nw.Physnetwork = Physnetwork
	#
	if visualize
		println("...  visualizing network structure ...")
		Visnetwork = Jumol.Network_vis(Dualnetwork, Physnetwork)
		nw.Visnetwork = Visnetwork
		Jumol.initialize_standard_plot_settings(nw.Visnetwork)
		Visnetwork.bfig = 2500
		Visnetwork.hfig = 2500
		plot_args = [true, true, false, true, true, false]
		Jumol.visualize_network(nw, plot_args)
	end
	# save in the struct
	lrs.nw = nw
end

## run LRS function
function run_lrs(lrs::LRS, x_coord_start=0.0, y_coord_start=0.0)
	println("... running LRS scan ...")
	## start the file
	output_fname = lrs.LI.output_path*"lrs_sig_"*string(lrs.SI.sig)*"_r_"*string(lrs.r_free)*"_sample_"*string(lrs.num_sample)*".dat"
	file_output = open(output_fname, "w")
	## generate the grid and prepare lists
	x_coord = x_coord_start
	y_coord = y_coord_start
	dx = lrs.nw.Dualnetwork.lx/lrs.LI.num_x_grid
	dy = lrs.nw.Dualnetwork.ly/lrs.LI.num_y_grid
	center_x = lrs.nw.Dualnetwork.lx*0.5
	center_y = lrs.nw.Dualnetwork.lx*0.5
	#
	lrs.x_coord_list = Vector()
	lrs.y_coord_list = Vector()
	for i in 1:lrs.LI.num_x_grid
		push!(lrs.x_coord_list, x_coord)
		x_coord += dx
	end
	for i in 1:lrs.LI.num_y_grid
		push!(lrs.y_coord_list, y_coord)
		y_coord += dy
	end
	#
	translate_vec = zeros(2)
	cntr_total = 1
	cntr_x = 1
	for x_coord in lrs.x_coord_list
		cntr_y = 1
		for y_coord in lrs.y_coord_list
			#
			println("cntr_x: ", cntr_x, ", cntr_y: ", cntr_y)
			write(file_output, "num_x_grid, num_y_grid\n")
			write(file_output, string(cntr_x)*" "*string(cntr_y)*"\n")
			write(file_output, "pos_x, pos_y\n")
			write(file_output, string(x_coord)*" "*string(y_coord)*"\n")
			# # translate the network to the center
			translate_vec[1] = center_x - x_coord
			translate_vec[2] = center_y - y_coord
			translate_network(lrs, translate_vec)
			dni, dbi, dti, pni, pbi = get_local_network_information(lrs)
			# save topological information and write into file
			mu, sig = extract_mean_ring_stat(lrs, dni)
			write(file_output, "ring statistics mu, sigma\n")
			write(file_output, string(mu)*" "*string(sig)*"\n")
			physical_nodes_mat = get_physical_nodes(lrs, pni)
			write(file_output, "physical nodes in matrix form: line1->x_coords, line2->y-coords\n")
			writedlm(file_output, physical_nodes_mat)
			bond_mat = calc_physical_bonds(lrs, pbi)
			write(file_output, "bond list in matrix form: line1->x_components, line2->y-components\n")
			writedlm(file_output, bond_mat)
			ring_Si_mat_list = calc_ring_Si_lists(lrs, dni)
			write(file_output, "Si ring nodes in matrix form: line1->x_coords, line2->y-coords (2 lines for one ring)\n")
			write(file_output, string(length(ring_Si_mat_list))*"\n")
			for ring_Si_mat in ring_Si_mat_list
				writedlm(file_output, ring_Si_mat)
			end
			# translate back to the original position
			translate_network(lrs, translate_vec*(-1.0))
			#
			cntr_y += 1
		end
		cntr_x += 1
	end
	close(file_output)
end

## build local dual and overlying local physical lattice
function get_local_network_information(lrs::LRS)
	# box center
	center_coord = zeros(2)
	center_coord[1] = lrs.nw.Dualnetwork.lx*0.5
	center_coord[2] = lrs.nw.Dualnetwork.ly*0.5
	# physical nodes
	phys_nodes_bool_list = Vector()
	for node in lrs.nw.Physnetwork.node_list
	    pos = zeros(2)
	    pos[1] = node.pos_x
	    pos[2] = node.pos_y
	    dist = norm(pos - center_coord)
	    if dist <= lrs.r_free
	        push!(phys_nodes_bool_list, true)
	    else
	        push!(phys_nodes_bool_list, false)
	    end
	end
	# physical bonds
	phys_bonds_bool_list = Vector()
	for bond in lrs.nw.Physnetwork.bond_list
	    if phys_nodes_bool_list[bond.node1.number] && phys_nodes_bool_list[bond.node2.number]
	        push!(phys_bonds_bool_list, true)
	    else
	        push!(phys_bonds_bool_list, false)
	    end
	end
	# dual triplets are equivalent to physical nodes
	dual_triplets_bool_list = phys_nodes_bool_list
	# dual nodes
	dual_node_bool_list = Vector()
	for dual_node in lrs.nw.Dualnetwork.dual_node_list
	    pos = zeros(2)
	    pos[1] = dual_node.pos_x
	    pos[2] = dual_node.pos_y
	    dist = norm(pos - center_coord)
	    if dist <= lrs.r_free
	        push!(dual_node_bool_list, true)
	    else
	        push!(dual_node_bool_list, false)
	    end
	end
	# dual bonds
	dual_bond_bool_list = Vector()
	for bond in lrs.nw.Dualnetwork.dual_bond_list
	    num1 = bond.node1.number
	    num2 = bond.node2.number
	    if dual_node_bool_list[num1] && dual_node_bool_list[num2]
	        push!(dual_bond_bool_list, true)
	    else
	        push!(dual_bond_bool_list, false)
	    end
	end
	return dual_node_bool_list, dual_bond_bool_list, dual_triplets_bool_list, phys_nodes_bool_list, phys_bonds_bool_list
end

## translate a whole network structure
function translate_network(lrs::LRS, translate_vec)
	# translate dual network
	for atom in lrs.nw.Dualnetwork.dualstruc.atom_list
	    atom.pos[1] += translate_vec[1]
	    atom.pos[2] += translate_vec[2]
	end
	put_atoms_back_to_box(lrs.nw.Dualnetwork.dualstruc)
	update_dual_network(lrs.nw.Dualnetwork)
	for triplet in lrs.nw.Dualnetwork.dual_triplet_list
	    find_triangle_nodes(triplet)
	end
	# translate physical network
	for atom in lrs.nw.Physnetwork.physstruc.atom_list
	    atom.pos[1] += translate_vec[1]
	    atom.pos[2] += translate_vec[2]
	end
	put_atoms_back_to_box(lrs.nw.Physnetwork.physstruc)
	update_node_pos(lrs.nw.Physnetwork)
	#translate silica network
	for atom in lrs.nw.Physnetwork.silicastruc.atom_list
	    atom.pos[1] += translate_vec[1]
	    atom.pos[2] += translate_vec[2]
	end
	put_atoms_back_to_box(lrs.nw.Physnetwork.silicastruc)
end

## plot local network structure
function plot_local_network(lrs::LRS, dni, dbi, dti, pni, pbi)
	println("...  visualizing local network structure ...")
	Visnetwork = Jumol.Network_vis(lrs.nw.Dualnetwork, lrs.nw.Physnetwork, false)
	lrs.nw.Visnetwork = Visnetwork
	Jumol.initialize_standard_plot_settings(Visnetwork)
	Visnetwork.bfig = 2500
	Visnetwork.hfig = 2500
	fig_local = Jumol.plot_box(Visnetwork)
	Jumol.plot_dual_nodes(Visnetwork, dni)
	Jumol.plot_dual_triplets(Visnetwork, dti)
	#Jumol.plot_bonds(Visnetwork, dbi)
	Jumol.plot_physical_bonds(Visnetwork, pbi)
	Jumol.plot_physical_nodes(Visnetwork, pni)
	Jumol.draw_circle(lrs, lrs.nw.Dualnetwork.lx*0.5, lrs.nw.Dualnetwork.ly*0.5, fig_local)
	display(fig_local)
	println("... press a button to continue ...")
	readline()
end

## draw the circle of the local region
function draw_circle(lrs::LRS, center_x, center_y, fig)
	N = 1000
	coords_x = zeros(N); coords_y = zeros(N); dphi = 2*pi/float(N)
	for i in 1:N
		coords_x[i] = center_x + lrs.r_free*cos(dphi*i)
		coords_y[i] = center_y + lrs.r_free*sin(dphi*i)
	end
	plot!(fig, coords_x, coords_y, color="magenta", linewidth=1.0, label=nothing)
end

## list of vectors to matrix
function vec2mat(vec_list)
	mat = zeros(length(vec_list[1]),length(vec_list))
	cntr = 1
	for vec in vec_list
		mat[:,cntr] = vec
		cntr += 1
	end
	return mat
end

## mean and variance of the ring statistics
function extract_mean_ring_stat(lrs, dni)
	## ring statistics matrix
	ring_count_mat = zeros(2,20)
	# fill first column (n-fold rings)
	for i = 1:size(ring_count_mat,2)
		ring_count_mat[1,i] = i
	end
	# fill second column (ring counts)
	cntr = 0
	for dn in lrs.nw.Dualnetwork.dual_node_list
		if dni[dn.number] == true
			ring_count_mat[2,dn.type] += 1
			cntr += 1
		end
	end
	# create ring stat matrix (area of the histogram is one)
	ring_stat_mat = zeros(2,20)
	for i in 1:size(ring_stat_mat,2)
		ring_stat_mat[1,i] = i
		ring_stat_mat[2,i] = float(ring_count_mat[2,i])/float(cntr)
	end
	## mean
	mu = 0.0
	for i in 1:size(ring_stat_mat,2)
		x = float(ring_stat_mat[1,i])
		y = float(ring_stat_mat[2,i])
		mu += y*x*1.0
	end
	mu = mu / 1.0 # area should be 1.0
	## standard deviation
	sig = 0.0
	for i in 1:size(ring_stat_mat,2)
		x = float(ring_stat_mat[1,i])
		y = float(ring_stat_mat[2,i])
		sig += (x-mu)*(x-mu)*y*1.0
	end
	sig = sqrt(sig/1.0)
	return mu, sig
end

## extract physical node list
function get_physical_nodes(lrs::LRS, pni)
	pnodes = Vector()
	cntr = 1
	for node in lrs.nw.Physnetwork.node_list
		if pni[cntr]
			coord = zeros(2)
			coord[1] = node.pos_x
			coord[2] = node.pos_y
			push!(pnodes, coord)
		end
		cntr += 1
	end
	pnodes_mat = vec2mat(pnodes)
	return pnodes_mat
end

## extract list of all bonds in the local region
function calc_physical_bonds(lrs::LRS, pbi)
	bond_list = Vector()
	cntr = 1
	for bond in lrs.nw.Physnetwork.bond_list
		if pbi[cntr]
			pos1 = [bond.node1.pos_x bond.node1.pos_y]
			pos2 = [bond.node2.pos_x bond.node2.pos_y]
			push!(bond_list, pos2-pos1)
		end
		cntr += 1
	end
	bond_mat = vec2mat(bond_list)
	return bond_mat
end

## extract Si atoms of every ring whose center is in the circle
function calc_ring_Si_lists(lrs::LRS, dni)
	cntr = 1
	Si_ring_list_list = Vector()
	for dn in lrs.nw.Dualnetwork.dual_node_list
		if dni[cntr]
			Si_num_list = dn.adjacent_triplet_list_ordered # because num_Si is equal to num_triplet
			Si_coord_list = Vector()
			for num in Si_num_list
				si_coord = [lrs.nw.Physnetwork.node_list[num].pos_x lrs.nw.Physnetwork.node_list[num].pos_y]
				push!(Si_coord_list, si_coord)
			end
			Si_coord_mat = vec2mat(Si_coord_list)
			push!(Si_ring_list_list, Si_coord_mat)
		end
		cntr += 1
	end
	return Si_ring_list_list
end

"""
function plot_ring_statistics(lrs::LRS)
	file_mu = lrs.LI.output_path*"mu_vis_sig"*string(lrs.SI.sig)*"_sample_"*string(lrs.num_sample)*".pdf"
	file_sig = lrs.LI.output_path*"sig_vis_sig"*string(lrs.SI.sig)*"_sample_"*string(lrs.num_sample)*".pdf"
	mu_rf = zeros(lrs.LI.num_x_grid, lrs.LI.num_y_grid)
	mu_bf = zeros(lrs.LI.num_x_grid, lrs.LI.num_y_grid)
	sig_rf = zeros(lrs.LI.num_x_grid, lrs.LI.num_y_grid)
	sig_bf = zeros(lrs.LI.num_x_grid, lrs.LI.num_y_grid)
	mu = zeros(lrs.LI.num_x_grid, lrs.LI.num_y_grid)
	sig = zeros(lrs.LI.num_x_grid, lrs.LI.num_y_grid)
	normalize_mu = zeros(lrs.LI.num_x_grid, lrs.LI.num_y_grid)
	normalize_sig = zeros(lrs.LI.num_x_grid, lrs.LI.num_y_grid)
	for i in 1:lrs.LI.num_x_grid
		for ii in 1:lrs.LI.num_y_grid
			ls = lrs.local_ring_stat[i][ii]
			mu[i,ii] = ls.local_mu
			sig[i,ii] = ls.local_sig
		end
	end
	max_mu = maximum(mu)
	min_mu = minimum(mu)
	max_sig = maximum(sig)
	min_sig = minimum(sig)
	for i in 1:lrs.LI.num_x_grid
		for ii in 1:lrs.LI.num_y_grid
			norm_mu = (mu[i,ii] - min_mu) / (max_mu - min_mu)
			norm_sig = (sig[i,ii] - min_sig) / (max_sig - min_sig)
			normalize_mu[i,ii] = norm_mu
			mu_bf[i,ii] = 1-norm_mu
			mu_rf[i,ii] = norm_mu
			normalize_sig[i,ii] = norm_sig
			sig_bf[i,ii] = 1-norm_sig
			sig_rf[i,ii] = norm_sig
		end
	end
	## plot and save mu
	fig_mu = plot(border=:none, aspect_ratio=1, legend=false, fmt=:pdf, title= "num_sample = "*string(lrs.num_sample)*", sig = "*string(lrs.SI.sig))
	plot!(size=(1000,1000))
	for i in 1:lrs.LI.num_x_grid
		for ii in 1:lrs.LI.num_y_grid
			ls = lrs.local_ring_stat[i][ii]
			x = ls.x_coord
			y = ls.y_coord
			fr = mu_rf[i,ii]
			fb = mu_bf[i,ii]
			scatter!([x],[y], markersize=13.0, markershape=:square,
	         		color=RGB(fr, 0.0, fb), markerstrokecolor=RGB(fr, 0.0, fb))
		end
	end
	plot_molecular_structure(lrs)
	savefig(fig_mu, file_mu)
	display(fig_mu)

	fig_sig = plot(border=:none, aspect_ratio=1, legend=false, fmt=:pdf, title= "num_sample = "*string(lrs.num_sample)*", sig = "*string(lrs.SI.sig))
	plot!(size=(lrs.PI.bfig,lrs.PI.hfig))
	for i in 1:lrs.LI.num_x_grid
		for ii in 1:lrs.LI.num_y_grid
			ls = lrs.local_ring_stat[i][ii]
			x = ls.x_coord
			y = ls.y_coord
			fr = sig_rf[i,ii]
			fb = sig_bf[i,ii]
			scatter!([x],[y], markersize=13.0, markershape=:square,
	         		color=RGB(fr, 0.0, fb), markerstrokecolor=RGB(fr, 0.0, fb))
		end
	end
	plot_molecular_structure(lrs)
	savefig(fig_sig, file_sig)
	display(fig_sig)
end

## plot the molecular structure
function plot_molecular_structure(lrs::LRS)
	molstruc = lrs.molstruc
	PI = lrs.PI
	Plotter = Vis2d(molstruc)
	Plotter.hfig=PI.hfig; Plotter.bfig=PI.bfig
	lx = molstruc.box.lx
	ly = molstruc.box.ly
	plot_atomic_structure_binary(Plotter, PI.at_size1, PI.at_size2,
	                             PI.factor_vec_x*lx,
								 PI.factor_vec_y*ly)
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
"""
