## Displacement field class
######################################
######################################
######################################
#
#
# (c) Franz Bamer Dec-2020
######################################
mutable struct Disp_field
    molstruc_init
    molstruc_step1
    molstruc_step2
    # limit value until the movement is plotted
    limit_val::Float64
    Disp_field(molstruc_init,
               molstruc_step1,
               molstruc_step2,
               limit_val = 0.1) = new(molstruc_init,
                                      molstruc_step1,
                                      molstruc_step2,
                                      limit_val)
    coord_mat_init::Array{Float64}
    coord_mat_step1::Array{Float64}
    coord_mat_step2::Array{Float64}
end


## calculate the coordinate matrices from the boxes
function get_coord_matrices(df::Disp_field)
    # init_box
    df.coord_mat_init = zeros(df.molstruc_init.noa,3)
    df.coord_mat_step1 = zeros(df.molstruc_step1.noa,3)
    df.coord_mat_step2 = zeros(df.molstruc_step2.noa,3)
    for i in 1:df.molstruc_init.noa
        df.coord_mat_init[i,:] = df.molstruc_init.atom_list[i].pos[1:3]
        df.coord_mat_step1[i,:] = df.molstruc_step1.atom_list[i].pos[1:3]
        df.coord_mat_step2[i,:] = df.molstruc_step2.atom_list[i].pos[1:3]
    end
end

## get affine displacement field from the boxes
function calc_aff_disp_field(df::Disp_field)
    ## triclinic box FIXME must be extended in z-direction
    aff_field = zeros(df.molstruc_init.noa,3)
    aff_field_init = zeros(df.molstruc_init.noa,3)
    #
    eps_x_1 = (df.molstruc_step1.box.lx - df.molstruc_init.box.lx) / df.molstruc_init.box.lx
    eps_x_2 = (df.molstruc_step2.box.lx - df.molstruc_init.box.lx) / df.molstruc_init.box.lx
    eps_x_12 = eps_x_2 - eps_x_1
    #
    eps_y_1 = (df.molstruc_step1.box.ly - df.molstruc_init.box.ly) / df.molstruc_init.box.ly
    eps_y_2 = (df.molstruc_step2.box.ly - df.molstruc_init.box.ly) / df.molstruc_init.box.ly
    eps_y_12 = eps_y_2 - eps_y_1
    #
    gamma_xy_0 = df.molstruc_init.box.lxy / df.molstruc_init.box.ly
    gamma_xy_1 = df.molstruc_step1.box.lxy / df.molstruc_step1.box.ly - gamma_xy_0
    gamma_xy_2 = df.molstruc_step2.box.lxy / df.molstruc_step2.box.ly - gamma_xy_0
    gamma_xy_12 = gamma_xy_2 - gamma_xy_1
    #
    for i in 1:df.molstruc_init.noa
        # affine displacement field x-direction
        aff_field[i,1] += eps_x_12 * df.molstruc_step1.atom_list[i].pos[1]
        # affine dislpacement field in y-direction
        aff_field[i,2] += eps_y_12 * df.molstruc_step1.atom_list[i].pos[2]
        # affine displacement field in xy-direction
        aff_field[i,1] += gamma_xy_12 * df.molstruc_step1.atom_list[i].pos[2]
        ##
        aff_field_init[i,1] = eps_x_1*df.molstruc_step1.atom_list[i].pos[1]
        aff_field_init[i,2] = eps_y_1*df.molstruc_step1.atom_list[i].pos[2]
        aff_field_init[i,1] += gamma_xy_1*df.molstruc_step1.atom_list[i].pos[2]
    end
    #
    return [aff_field, aff_field_init]
end

## get affine displacement field from the boxes
function calc_aff_disp_field_lmp(df::Disp_field)
    ## triclinic box FIXME must be extended in z-direction
    aff_field = zeros(df.molstruc_init.noa,3)
    aff_field_init = zeros(df.molstruc_init.noa,3)
    #
    eps_x_12 = (df.molstruc_step2.box.lx - df.molstruc_step1.box.lx) / df.molstruc_step1.box.lx
    #
    eps_y_12 = (df.molstruc_step2.box.ly - df.molstruc_step1.box.ly) / df.molstruc_step1.box.ly
    #
    gamma_xy_0 = df.molstruc_init.box.lxy / df.molstruc_init.box.ly
    gamma_xy_1 = df.molstruc_step1.box.lxy / df.molstruc_step1.box.ly - gamma_xy_0
    gamma_xy_2 = df.molstruc_step2.box.lxy / df.molstruc_step2.box.ly - gamma_xy_0
    gamma_xy_12 = gamma_xy_2 - gamma_xy_1
    #
    for i in 1:df.molstruc_init.noa
        # affine displacement field x-direction
        aff_field[i,1] += eps_x_12 * df.molstruc_step1.atom_list[i].pos[1]
        # affine dislpacement field in y-direction
        aff_field[i,2] += eps_y_12 * df.molstruc_step1.atom_list[i].pos[2]
        # affine displacement field in xy-direction
        aff_field[i,1] += gamma_xy_12 * df.molstruc_step1.atom_list[i].pos[2]
        ##
        #aff_field_init[i,1] = eps_x_1*df.molstruc_step1.atom_list[i].pos[1]
        #aff_field_init[i,2] = eps_y_1*df.molstruc_step1.atom_list[i].pos[2]
        #aff_field_init[i,1] += gamma_xy_1*df.molstruc_step1.atom_list[i].pos[2]
    end
    #
    return [aff_field, aff_field_init]
end

## get full displacement field
function calc_full_disp_field(df::Disp_field)
    field_total = df.coord_mat_step2 - df.coord_mat_step1
    return field_total
end

## get non-affine displacement field
function get_non_aff_disp_field(df::Disp_field)

end

## plot a dislpacement field
function plot_displacement_field(df::Disp_field, field::Array{Float64},
                                 field_init::Array{Float64}, factor::Float64=1.0, factor_lw::Float64=10.0)
    coord_mat_init = df.coord_mat_step1# - field_init #FIXME only works for simple shear
    # quiver plot of the input field
    for i=1:df.molstruc_init.noa
        #sze = norm([field[i,1],field[i,2]])*factor
        #quiver!([coord_mat_init[i,1]], [coord_mat_init[i,2]],
        #        quiver=([field[i,1]].*factor, [field[i,2]].*factor), linewidths=sze,
        #        color="black")
        plot_arrow(df, coord_mat_init[i,1:2], field[i,1:2], factor, factor_lw)
    end
end


## custom arrow function plot
function plot_arrow(df::Disp_field, start_coord::Vector,
                    vector::Vector, factor::Float64, factor_lw::Float64)
    # displacement vector
    len = norm(vector)
    if len < df.limit_val
        # resize length of vector
        len = len*factor
        linew = len/df.molstruc_init.box.lx * factor_lw
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
