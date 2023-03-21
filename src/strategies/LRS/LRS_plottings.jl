## Read LRS results from file
#  plot functions
######################################
######################################
######################################
#
#
# (c) Franz Bamer, Ivy Wu September-2022
######################################






## plot local structure
function plot_local_structure_rings(lr::LRS_Read, r_free, index, center_x, center_y)
	println(" ... plotting local structure ...")
	println("pos_x: ", lr.pos_list[index][1], " pos_y: ", lr.pos_list[index][2])
	#
	fig_local = plot(size=(1000,1000), aspect_ratio=1,
	                 title="local structure; pos_x="*string(lr.pos_list[index][1])*"; pos_y="*string(lr.pos_list[index][2]), label=nothing)
	#
	for ring_Si_list in lr.ring_Si_list_list_list[index]
		coords_x = Vector(); coords_y = Vector()
		for Si_atom in ring_Si_list
			push!(coords_x, Si_atom[1])
			push!(coords_y, Si_atom[2])
		end
		push!(coords_x, ring_Si_list[1][1])
		push!(coords_y, ring_Si_list[1][2])
		# plot ring
		plot!(fig_local, coords_x, coords_y, linewidth=1.0, color="blue", label=nothing)
		scatter!(fig_local, coords_x, coords_y, markersize=4.0, color="blue", label=nothing)
	end
	draw_circle(lr, r_free, center_x, center_y, fig_local)
	return fig_local
end

## draw the circle of the local region
function draw_circle(lr::LRS_Read, r_free, center_x, center_y, fig)
	N = 1000
	coords_x = zeros(N); coords_y = zeros(N); dphi = 2*pi/float(N)
	for i in 1:N
		coords_x[i] = center_x + r_free*cos(dphi*i)
		coords_y[i] = center_y + r_free*sin(dphi*i)
	end
	plot!(fig, coords_x, coords_y, color="magenta", linewidth=1.0, label=nothing)
end

## customized arrow function plot
function plot_arrow(start_coord, vector, plot_factor, text=" ", col="blue")
	#
	vector = vector * plot_factor
	#
	end_coord = start_coord + vector
    # displacement vector
    len = norm(vector)
    # rotation matrix
    cos_phi = vector[1]/len
    sin_phi = vector[2]/len
    T = [cos_phi -sin_phi
         sin_phi  cos_phi]
    # arrow head
    coord_arrow_head = [ -len*0.25  0.0  -len*0.25
                          len*0.15  0.0  -len*0.15 ]
    # rotate arrow head
    coord_arrow_head = T*coord_arrow_head
    # translate arrow head
    coord_arrow_head[1,:] = coord_arrow_head[1,:] .+ end_coord[1]
    coord_arrow_head[2,:] = coord_arrow_head[2,:] .+ end_coord[2]
    # plot arrow
    plot!([start_coord[1],end_coord[1]], [start_coord[2],end_coord[2]],
          linewidth=3.0, color=col, label=nothing)
    # plot arrow head
    plot!(coord_arrow_head[1,:],coord_arrow_head[2,:],
          linewidth=3.0, color=col, label=nothing)
	annotate!(coord_arrow_head[1,2],coord_arrow_head[2,2], text, color=col)
end


function plot_bond_vector_fig(lr::LRS_Read, index::Int64)
	println("... plotting bond vectors ...")
	#
	fig_bonds = plot(size=(1000,1000), aspect_ratio=1,
	                 title="bond vectors", label=nothing)
	bond_list = mat2vec_list(lr.bond_mat_list[index])
	for bond in bond_list
		plot!(fig_bonds, [0,bond[1]],[0,bond[2]], color="black", linewidth=2.0, label=nothing)
	end
	return fig_bonds
end

function plot_radial_hist(lr::LRS_Read, hist)
	dphi = 2.0*pi/float(lr.num_intersec_hist)
	# draw historgram
	N = 200
	dphi_s = dphi/float(N)
	fig_hist = plot(size=(1000,1000), aspect_ratio=1, xlims=(-0.5,0.5), ylims=(-0.5,0.5),
	                 title="histogram", label=nothing)
	plot!([-1,1],[0,0], color="black", linewidth=1.0, linestyle=:dot, label=nothing)
	plot!([0,0],[-1,1], color="black", linewidth=1.0, linestyle=:dot, label=nothing)
	coords_x = zeros(N*lr.num_intersec_hist+1)
	coords_y = zeros(N*lr.num_intersec_hist+1)
	cntr = 1
	for i in 1:lr.num_intersec_hist
		angle_start = (i-1) * dphi
		for ii in 1:N
			coords_x[cntr] = hist[i]*cos(angle_start + dphi_s*(ii-1))
			coords_y[cntr] = hist[i]*sin(angle_start + dphi_s*(ii-1))
			cntr += 1
		end
	end
	coords_x[cntr] = coords_x[1]
	coords_y[cntr] = coords_y[1]
	plot!(fig_hist, coords_x, coords_y, color="magenta", linewidth=2.0, label=nothing)
	#
	#for bond in bond_list
	#	plot!(fig_hist, [0,bond[1]],[0,bond[2]], color="black", linewidth=2.0, label=nothing)
	#end
	#
	return fig_hist
end
