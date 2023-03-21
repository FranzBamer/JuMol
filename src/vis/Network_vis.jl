## plot the dual structure (2d)
######################################
######################################
######################################




using Plots
using PrettyTables


mutable struct Network_vis
    dns # dual network structure
	pns # physical network structure
	pbc
    Network_vis(dns, pns, pbc=true) = new(dns, pns, pbc)
	#
	dual_node_loc_info::Vector
	dual_bond_loc_info::Vector
	dual_triplet_loc_infor::Vector
	node_loc_info::Vector
	bond_loc_info::Vector
	## figure size settings
	bfig::Float64
	hfig::Float64
	## box settings
	box_color::String
	box_lw::Float64
	## dual node settings
	dual_node_color::String
	dual_node_numbering::Bool
	dual_node_size::Float64
	dual_node_text_size::Int64
	dual_node_plot_x_shift::Float64
	dual_node_plot_y_shift::Float64
	## dual triplet settings
	dual_triplet_color::String
	dual_triplet_numbering::Bool
	dual_triplet_text_size::Int64
	dual_triplet_lw::Float64
	## dual bond settings
	dual_bond_color::String
	dual_bond_textsize::Int64
	dual_bond_numbering::Bool
	dual_bond_lw::Float64
	dual_bond_text_x_shift::Float64
	dual_bond_text_y_shift::Float64
	## physical network nodes settings
	node_color::String
	node_numbering::Bool
	node_size::Float64
	node_textsize::Int64
	node_text_x_shift::Float64
	node_text_y_shift::Float64
	## physical bond settings
	bond_color::String
	bond_lw::Float64
	bond_numbering::Bool
	bond_textsize::Int64
	bond_text_x_shift::Float64
	bond_text_y_shift::Float64
	## physical bond node settings
	node_bond_color::String
	node_bond_size::Float64
	## silica bond settings
	silica_bond_color::String
	silica_bond_numbering::Bool
	silica_text_size::Int64
	silica_bond_lw::Float64
	silica_bond_text_x_shift::Float64
	silica_bond_text_y_shift::Float64
end

## initialize plot settings (always call this before plotting)
function initialize_standard_plot_settings(nv::Network_vis)
	## figure size
	nv.bfig = 500.0
	nv.hfig = 500.0
	## box settings
	nv.box_color = "gray"
	nv.box_lw = 2.0
	## dual node settings
	nv.dual_node_color = "red"
	nv.dual_node_numbering = false
	nv.dual_node_size = 4
	nv.dual_node_text_size = 8
	nv.dual_node_plot_x_shift = -0.3
	nv.dual_node_plot_y_shift = -0.1
	## dual triplet settings
	nv.dual_triplet_color = "black"
	nv.dual_triplet_numbering = false
	nv.dual_triplet_text_size = 6
	nv.dual_triplet_lw = 2.0
	## dual bond settings
	nv.dual_bond_color = "gray"
	nv.dual_bond_lw = 1.0
	nv.dual_bond_textsize = 10
	nv.dual_bond_numbering = false
	nv.dual_bond_text_x_shift = -0.3
	nv.dual_bond_text_y_shift = -0.1
	## physical node settings
	nv.node_color = "blue"
	nv.node_numbering = false
	nv.node_size = 5.0
	nv.node_textsize = 8
	nv.node_text_x_shift = -0.3
	nv.node_text_y_shift = -0.1
	## physical bond settings
	nv.bond_color = "blue"
	nv.bond_lw = 2.0
	nv.bond_numbering = false
	nv.bond_textsize = 10
	nv.bond_text_x_shift = 0.0
	nv.bond_text_y_shift = 0.0
	## physical bond node settings
	nv.node_bond_color = "blue"
	nv.node_bond_size = 0.0
	## silica bond settings
	nv.silica_bond_color = "blue"
	nv.silica_bond_numbering = false
	nv.silica_text_size = 10
	nv.silica_bond_lw = 2.0
	nv.silica_bond_text_x_shift = 0.0
	nv.silica_bond_text_y_shift = 0.0
end


## plot the box
function plot_box(nv::Network_vis)
	p_mat = [ 0                     0
	          nv.dns.lx  0
			  nv.dns.lx  nv.dns.ly
			  0                     nv.dns.ly
			  0                     0                       ]
	box_border_plot = plot(p_mat[:,1],p_mat[:,2],border=:none,aspect_ratio=1,
	            legend=false,color=nv.box_color,linewidth=nv.box_lw,fmt=:pdf)
	plot!(size=(nv.bfig,nv.hfig))
	return box_border_plot
end

## plot dual nodes
function plot_dual_nodes(nv::Network_vis, loc_info=Vector())
	coords_x = Vector()
	coords_y = Vector()
	if length(loc_info) == length(nv.dns.dual_node_list)
		cntr = 1
		for n in nv.dns.dual_node_list
			if loc_info[cntr]
				push!(coords_x, n.pos_x)
				push!(coords_y, n.pos_y)
			end
			cntr += 1
		end
	else
		for n in nv.dns.dual_node_list
			push!(coords_x, n.pos_x)
			push!(coords_y, n.pos_y)
		end
	end
	scatter!(coords_x, coords_y,
		     markersize=nv.dual_node_size,color=nv.dual_node_color,
			 markerstrokecolor=nv.dual_node_color)
	if nv.dual_node_numbering == true
		for i in 1:num_nodes
			n = nv.dns.dual_node_list[i]
			annotate!(n.pos_x+nv.dual_node_plot_x_shift,
			          n.pos_y+nv.dual_node_plot_y_shift,
			          text(string(i), :center,
					  color=nv.dual_node_color, nv.dual_node_text_size))
		end
	end
end

## plot dual triplets
function plot_dual_triplets(nv::Network_vis, loc_info=Vector())
	#
	for i in 1:length(nv.dns.dual_triplet_list)
		tr = nv.dns.dual_triplet_list[i]
		num1 = tr.node1.number
		num2 = tr.node2.number
		num3 = tr.node3.number
		if nv.pbc
			pos1 = nv.dns.dualstruc.atom_list[num1].pos[1:2]
			pos2 = pos1 - get_distance_between2atoms(nv.dns.dualstruc,num1,num2)[1:2]
			pos3 = pos2 - get_distance_between2atoms(nv.dns.dualstruc,num2,num3)[1:2]
		else
			pos1 = [nv.dns.dual_node_list[num1].pos_x nv.dns.dual_node_list[num1].pos_y]
			pos2 = [nv.dns.dual_node_list[num2].pos_x nv.dns.dual_node_list[num2].pos_y]
			pos3 = [nv.dns.dual_node_list[num3].pos_x nv.dns.dual_node_list[num3].pos_y]
		end
		##
		triplet_mat = [ pos1[1] pos2[1] pos3[1] pos1[1]
		                pos1[2] pos2[2] pos3[2] pos1[2] ]
		if length(loc_info) == length(nv.dns.dual_triplet_list)
			if loc_info[i]
				plot!(triplet_mat[1,:],triplet_mat[2,:],
				      color=nv.dual_triplet_color,linewidth=nv.dual_triplet_lw)
				if nv.dual_triplet_numbering==true
					pos_center = (pos1 + pos2 + pos3)/3.0
					annotate!(pos_center[1],pos_center[2],
					          text(string(i), color=nv.dual_triplet_color, :center,
							  nv.dual_triplet_text_size))
				end
			end
		else
			plot!(triplet_mat[1,:],triplet_mat[2,:],
			      color=nv.dual_triplet_color,linewidth=nv.dual_triplet_lw)
			if nv.dual_triplet_numbering==true
				pos_center = (pos1 + pos2 + pos3)/3.0
				annotate!(pos_center[1],pos_center[2],
				          text(string(i), color=nv.dual_triplet_color, :center,
						  nv.dual_triplet_text_size))
			end
		end
	end
end

## plot dual bonds
function plot_bonds(nv::Network_vis, loc_info=Vector())
	for i in 1:length(nv.dns.dual_bond_list)
		b = nv.dns.dual_bond_list[i]
		num1 = b.node1.number
		num2 = b.node2.number
		if nv.pbc
			pos1 = nv.dns.dualstruc.atom_list[num1].pos[1:2]
			pos2 = pos1 - get_distance_between2atoms(nv.dns.dualstruc,num1,num2)[1:2]
		else
			pos1 = [nv.dns.dual_node_list[num1].pos_x nv.dns.dual_node_list[num1].pos_y]
			pos2 = [nv.dns.dual_node_list[num2].pos_x nv.dns.dual_node_list[num2].pos_y]
		end
		bond_mat = [ pos1[1] pos2[1]
		             pos1[2] pos2[2]]
		if length(loc_info) == length(nv.dns.dual_bond_list)
			if loc_info[i]
				plot!(bond_mat[1,:],bond_mat[2,:],
				      color=nv.dual_bond_color,linewidth=nv.dual_bond_lw)
				if nv.dual_bond_numbering == true
					pos_center_bond = pos1 + 0.5*(pos2-pos1)
					annotate!(pos_center_bond[1]+nv.dual_bond_text_x_shift,
					          pos_center_bond[2]+nv.dual_bond_text_y_shift,
							  text(string(i), color=nv.dual_bond_color, :center,
							  nv.dual_bond_textsize))
				end
			end
		else
			plot!(bond_mat[1,:],bond_mat[2,:],
			      color=nv.dual_bond_color,linewidth=nv.dual_bond_lw)
			if nv.dual_bond_numbering == true
				pos_center_bond = pos1 + 0.5*(pos2-pos1)
				annotate!(pos_center_bond[1]+nv.dual_bond_text_x_shift,
				          pos_center_bond[2]+nv.dual_bond_text_y_shift,
						  text(string(i), color=nv.dual_bond_color, :center,
						  nv.dual_bond_textsize))
			end
		end
	end
end

## plot network nodes
function plot_physical_nodes(nv::Network_vis, loc_info=Vector())
	coords_x = Vector()
	coords_y = Vector()
	if length(loc_info) == length(nv.pns.node_list)
		cntr = 1
		for n in nv.pns.node_list
			if loc_info[cntr]
				push!(coords_x, n.pos_x)
				push!(coords_y, n.pos_y)
			end
			cntr += 1
		end
	else
		for i in 1:length(nv.pns.node_list)
			n = nv.pns.node_list[i]
			push!(coords_x, n.pos_x)
			push!(coords_y, n.pos_y)
		end
	end
	scatter!(coords_x,coords_y,
		     markersize=nv.node_size,color=nv.node_color,
			 markerstrokecolor=nv.node_color)
	if nv.node_numbering == true
		for i in 1:num_nodes
			n = nv.pns.node_list[i]
			annotate!(n.pos_x+nv.node_text_x_shift,n.pos_y+nv.node_text_y_shift,
			          text(string(i), color=nv.node_color, :center,
					  nv.node_textsize))
		end
	end
end

## plot network bonds
function plot_physical_bonds(nv::Network_vis, loc_info=Vector())
	for i in 1:length(nv.pns.bond_list)
		b = nv.pns.bond_list[i]
		num1 = b.node1.number
		num2 = b.node2.number
		if nv.pbc
			pos1 = nv.pns.physstruc.atom_list[num1].pos[1:2]
			pos2 = pos1 - get_distance_between2atoms(nv.pns.physstruc,num1,num2)[1:2]
		else
			pos1 = [nv.pns.node_list[num1].pos_x nv.pns.node_list[num1].pos_y]
			pos2 = [nv.pns.node_list[num2].pos_x nv.pns.node_list[num2].pos_y]
		end
		bond_mat = [ pos1[1] pos2[1]
		             pos1[2] pos2[2] ]
		if length(loc_info) == length(nv.pns.bond_list)
			if loc_info[i]
				plot!(bond_mat[1,:],bond_mat[2,:],
				      color=nv.bond_color,linewidth=nv.bond_lw)
				if nv.bond_numbering == true
					pos_center_bond = pos1 + 0.5*(pos2-pos1)
					annotate!(pos_center_bond[1]+nv.bond_text_x_shift,
					          pos_center_bond[2]+nv.bond_text_y_shift,
							  text(string(i), color=nv.bond_color, :center,
							  nv.bond_textsize))
				end
			end
		else
			plot!(bond_mat[1,:],bond_mat[2,:],
			      color=nv.bond_color,linewidth=nv.bond_lw)
			if nv.bond_numbering == true
				pos_center_bond = pos1 + 0.5*(pos2-pos1)
				annotate!(pos_center_bond[1]+nv.bond_text_x_shift,
				          pos_center_bond[2]+nv.bond_text_y_shift,
						  text(string(i), color=nv.bond_color, :center,
						  nv.bond_textsize))
			end
		end
	end
end

## plot silica network bonds
function plot_physical_silica_bonds(nv::Network_vis)
	for i in 1:length(nv.pns.bond_list)
		b = nv.pns.bond_list[i]
		num1 = b.node1.number
		num2 = b.node2.number
		pos1 = nv.pns.silicastruc.atom_list[num1].pos[1:2]
		pos2 = pos1 - get_distance_between2atoms(nv.pns.silicastruc,num1,num2)[1:2]
		bond_mat = [ pos1[1] pos2[1]
		             pos1[2] pos2[2] ]
		# plot the bond
		plot!(bond_mat[1,:],bond_mat[2,:],
		      color=nv.silica_bond_color,linewidth=nv.silica_bond_lw)
		# plot the bond atoms
		pos_center_bond = pos1 + 0.5*(pos2-pos1)
		scatter!([pos_center_bond[1]],[pos_center_bond[2]],
			     markersize=nv.node_bond_size,color=nv.node_bond_color,
				 markerstrokecolor=nv.node_bond_color)
		if nv.silica_bond_numbering == true
			annotate!(pos_center_bond[1]+nv.silica_bond_text_x_shift,
			          pos_center_bond[2]+nv.silica_bond_text_y_shift,
					  text(string(i), color=nv.silica_bond_color, :center,
					  nv.silica_bond_textsize))
		end
	end
end
