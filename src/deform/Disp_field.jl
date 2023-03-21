## Displacement field class
######################################
######################################
######################################
#
#
# (c) Franz Bamer Dec-2020
######################################

using LinearAlgebra

mutable struct Disp_field
    molstruc_step1
    molstruc_step2
    # limit value until the movement is plotted
    limit_val::Float64
    Disp_field(molstruc_step1,
               molstruc_step2,
               limit_val = 0.1) = new(molstruc_step1,
                                      molstruc_step2,
                                      limit_val)
    coord_mat_step1::Array{Float64}
    coord_mat_step2::Array{Float64}
    coord_mat_unit_square::Array{Float64}
    #
    H_1::Array{Float64}
    H_2::Array{Float64}
	#
	aff_field::Array{Float64}
	field_total::Array{Float64}
	nonaff_field::Array{Float64}
end

## get affine displacement field from the boxes
function calc_disp_fields(df::Disp_field)
	# coordinate matrices
	df.coord_mat_step1 = zeros(df.molstruc_step1.noa,3)
    df.coord_mat_step2 = zeros(df.molstruc_step2.noa,3)
    for i in 1:df.molstruc_step1.noa
        df.coord_mat_step1[i,:] = df.molstruc_step1.atom_list[i].pos[1:3]
        df.coord_mat_step2[i,:] = df.molstruc_step2.atom_list[i].pos[1:3]
    end
    # triclinic box at step 1
    df.H_1 = zeros(3,3)
    df.H_1[:,1] = df.molstruc_step1.box.h1
    df.H_1[:,2] = df.molstruc_step1.box.h2
    df.H_1[:,3] = df.molstruc_step1.box.h3
	#println("H_1")
	#println(df.H_1)
    # triclinic box at step 2
    df.H_2 = zeros(3,3)
    df.H_2[:,1] = df.molstruc_step2.box.h1
    df.H_2[:,2] = df.molstruc_step2.box.h2
    df.H_2[:,3] = df.molstruc_step2.box.h3
	#println("H_2")
	#println(df.H_2)
	# map the coordinates from step 1 to the unit box
	df.coord_mat_unit_square = zeros(df.molstruc_step1.noa,3)
	H_1_inv = inv(df.H_1)
	#println("H_1_inv")
	#println(H_1_inv)
    for i in 1:df.molstruc_step1.noa
        df.coord_mat_unit_square[i,:] = H_1_inv*df.coord_mat_step1[i,:]
    end
    # deformation gradient and displacement gradient from step1 to step2
    F_12 = df.H_2*inv(df.H_1)
	#println("F_12")
	#println(F_12)
	D_12 = F_12 - Matrix{Float64}(I,3,3)
	#println("D_12")
	#println(D_12)
    # affine displacement field for every atom
    df.aff_field = zeros(df.molstruc_step1.noa,3)
    for i in 1:df.molstruc_step1.noa
        df.aff_field[i,:] = D_12*df.molstruc_step1.atom_list[i].pos
    end
	# full displacement field
	df.field_total = df.coord_mat_step2 - df.coord_mat_step1
	# non-affine displacement field
	df.nonaff_field = df.field_total - df.aff_field
	#println("aff field")
	#println(df.aff_field[1:10,:])
	#println("field total")
	#println(df.field_total[1:10,:])
	#println("non-aff field")
	#println(df.nonaff_field[1:10,:])
end

## get atom index and coordinate that experiences the maximum displacement between the steps
function locate_max_displacement(df::Disp_field, maximum_limit_disp)
	max_disp = -1.0
	atom_index = 0
	#println("disp:")
	for i in 1:size(df.nonaff_field,1)
		atomic_disp_vec = df.nonaff_field[i,:]
		disp = norm(atomic_disp_vec)
		#if i < 10
		#	print(disp,", ")
		#end
		if disp > max_disp && disp < maximum_limit_disp
			max_disp = disp
			atom_index = i
		end
	end
	#println("maximum found")
	#println(atom_index)
	#readline()
	return atom_index
end

## plot the unit box
function plot_unit_box(df::Disp_field, linewidth=0.5)
    p_mat = zeros(3,5)
	p_mat[:,2] = zeros(3) + [1;0;0]
	p_mat[:,3] = p_mat[:,2] + [0;1;0]
	p_mat[:,4] = p_mat[:,3] - [1;0;0]
	unit_box_plot = plot(p_mat[1,:],p_mat[2,:],border=:none,aspect_ratio=1,
	            legend=false,color="black",lw=linewidth,fmt=:pdf)
	plot!(size=(500,500))
	return unit_box_plot
end

## plot a displacement field
# create a plotting object before calling the function
function plot_displacement_field(df::Disp_field, field::Array{Float64}, boxlw::Float64,
                                 factor::Float64=1.0, factor_lw::Float64=10.0)
	p_mat = zeros(3,5)
	# plot unit box borders
	p_mat[:,2] = zeros(3) + [1;0;0]
	p_mat[:,3] = p_mat[:,2] + [0;1;0]
	p_mat[:,4] = p_mat[:,3] - [1;0;0]
	fig = plot(p_mat[1,:],p_mat[2,:],border=:none,aspect_ratio=1,
	            legend=false,color="black",lw=boxlw, fmt=:pdf)
	plot!(size=(500,500))
	factor = factor * 1.0/(df.molstruc_step1.box.lx)
	# quiver plot of the input field
    for i=1:df.molstruc_step1.noa
        plot_arrow(df, df.coord_mat_unit_square[i,1:2], field[i,1:2], factor, factor_lw)
    end
	display(fig)
	return fig
end

## custom arrow function plot
function plot_arrow(df::Disp_field, start_coord::Vector,
                    vector::Vector, factor::Float64, factor_lw::Float64)
    # displacement vector
    len = norm(vector)
    if len < df.limit_val
        # resize length of vector
        len = len*factor
        linew = len * factor_lw
        vector = vector .* factor
        end_coord = start_coord + vector
        # rotation matrix
        cos_phi = vector[1]/len
        sin_phi = vector[2]/len
        T = [cos_phi -sin_phi
             sin_phi  cos_phi]
        # arrow head
        coord_arrow_head = [ -len*0.35  0.0  -len*0.35
                              len*0.2  0.0  -len*0.2 ]
        # rotate arrow head
        coord_arrow_head = T*coord_arrow_head
        # translate arrow head
        coord_arrow_head[1,:] = coord_arrow_head[1,:] .+ end_coord[1]
        coord_arrow_head[2,:] = coord_arrow_head[2,:] .+ end_coord[2]

        # plot arrow
        plot!([start_coord[1],end_coord[1]], [start_coord[2],end_coord[2]],
              linewidth=linew, color="black")
        # plot arrow head
        plot!(coord_arrow_head[1,:],coord_arrow_head[2,:],
              linewidth=linew, color="black")
    end
end
