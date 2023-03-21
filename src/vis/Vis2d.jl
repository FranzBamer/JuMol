## plot structure (2d)
######################################
######################################
######################################




using Plots
using PrettyTables





mutable struct Vis2d
	structure
	bfig::Float64
	hfig::Float64
	box_border_plot
	Vis2d(structure,bfig=500,hfig=500) = new(structure,bfig,hfig)
end

## must be called before every plot
function plot_box(plotter::Vis2d, box_color="gray",linewidth=1.0)
	p_mat = zeros(3,5)
	p_mat[:,2] = zeros(3) + plotter.structure.box.h1
	p_mat[:,3] = p_mat[:,2] + plotter.structure.box.h2
	p_mat[:,4] = p_mat[:,3] - plotter.structure.box.h1
	box_border_plot = plot(p_mat[1,:],p_mat[2,:],border=:none,aspect_ratio=1,
	            legend=false,color=box_color,lw=linewidth,fmt=:pdf)
	plot!(size=(plotter.bfig,plotter.hfig))
	return box_border_plot
end

## can be called at the end of the plot if wanted
function plot_box_border(plotter::Vis2d, box_color, linewidth)
	p_mat = zeros(3,5)
	p_mat[:,2] = zeros(3) + plotter.structure.box.h1
	p_mat[:,3] = p_mat[:,2] + plotter.structure.box.h2
	p_mat[:,4] = p_mat[:,3] - plotter.structure.box.h1
	box_border_plot = plot!(p_mat[1,:],p_mat[2,:],
	                        color=box_color, lw=linewidth)
	plot!(size=(plotter.bfig,plotter.hfig))
end

function plot_atomic_structure_binary(plotter::Vis2d, size1::Float64, size2::Float64,
                               framex::Vector, framey::Vector, col1="lightblue", col2="lightgray",
							   linewidth=0.5, op=1.0)
    Jumol.plot_white_frame(plotter,framex, framey)
    Jumol.plot_atoms_type(plotter, 1, "black", col1, size1, linewidth, op)
    Jumol.plot_atoms_type(plotter, 2, "black", col2, size2, linewidth, op)
end

function plot_covalent_bonds(plotter::Vis2d, bond_color="blue", bondwidth=1, cov_dist=1.8)
	for i in 1:plotter.structure.noa
		atom = plotter.structure.atom_list[i]
		for ii in atom.neighbor_indices
			if ii > i
				dist = get_distance_between2atoms(plotter.structure, i, ii)
				if dist[4] < cov_dist
					plot_mat = zeros(2,2)
					plot_mat[:,1] = atom.pos[1:2]
					plot_mat[:,2] = atom.pos[1:2] - dist[1:2]
					plot!(plot_mat[1,:],plot_mat[2,:],
					      linewidth=bondwidth, color=bond_color)
				end
			end
		end
	end
end

function plot_atoms(plotter::Vis2d,atom_color="blue",atom_size=10)
	scatter!(get_atom_pos_mat(plotter.structure)[:,1],
	      get_atom_pos_mat(plotter.structure)[:,2],
		  markersize=atom_size,color=atom_color,markerstrokecolor=atom_color)
end

function plot_atoms_type(plotter::Vis2d, atom_type::Int64, atom_color="blue",
	                     fill_color="white", atom_size=10, linewidth=1.0, op=1.0)
	for atom in plotter.structure.atom_list
		if atom.type == atom_type
			if atom.group == 0
				scatter!([atom.pos[1]], [atom.pos[2]],
				         markersize=atom_size, color=fill_color, markerstrokecolor=atom_color, markerstrokewidth=linewidth, opacity=op)
			else
				scatter!([atom.pos[1]], [atom.pos[2]],
				         markersize=atom_size, color="gray", markerstrokecolor="gray", markerstrokewidth=linewidth, opacity=op)
			end
		end
	end
end

function plot_white_frame(plotter::Vis2d,x_coords,y_coords)
	plot!(x_coords,y_coords,color="white",lw=1)
end

function plot_coordinate(plotter::Vis2d, size, x_coord, y_coord)
	scatter!([x_coord], [y_coord],
			 marker=:x, markersize=size, color="magenta", markerstrokecolor="magenta", markerstrokewidth=0.1)
end

function plot_circle(plotter::Vis2d, center_x, center_y, radius)
	N = 1000
	x_coords = zeros(N)
	y_coords = zeros(N)
	dphi = 2.0*pi / float(N)
	for i in 1:N
		x_coords[i] = center_x + radius * cos(dphi*i)
		y_coords[i] = center_y + radius * sin(dphi*i)
	end
	plot!(x_coords, y_coords, color="black", linewidth=2.0)
end

function build_tikz_picture(plotter::Vis2d, fname="atom_struc.tex")
	f = open(fname, "w")
	write(f,"starting line\n")
	close(f)
end
