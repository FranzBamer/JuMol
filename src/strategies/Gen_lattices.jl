## create lattice structures
######################################
######################################
######################################



mutable struct Gen_lattice
	molstruc
	Gen_lattice(molstruc) = new(molstruc)
end

## create a LJ lattice structure -- hexagonally packed structure
function create_LJ_lattice(gl::Gen_lattice, num_x::Int64, num_y::Int64)
        c = 0.693719696937019
        create_box_by_hand(gl.molstruc,num_x*c,(num_y*2)*c*sin(pi/3),10.0,0.0,0.0,0.0)
        #
        atom_cntr = 1
        #
        for i in 1:num_y
                x_start1 = 0.0
                y_start1 = 2*(i-1)*c*sin(pi/3)
                for ii in 1:num_x
                        Jumol.add_atom_by_hand(gl.molstruc,atom_cntr,1,x_start1+(ii-1)*c,y_start1,0.0,0.0,0.0,0.0)
                        atom_cntr += 1
                end
                x_start2 = c*cos(pi/3)
                y_start2 = y_start1 + c*sin(pi/3)
                for ii in 1:num_x
                        Jumol.add_atom_by_hand(gl.molstruc,atom_cntr,1,x_start2+(ii-1)*c,y_start2,0.0,0.0,0.0,0.0)
                        atom_cntr += 1
                end
        end
end

## generate a LJ glass configuration ratio small to large is the golden mean
function create_random_pos_LJ(gl::Gen_lattice, num_x, num_y)
        lbox_unit = 0.988045
        lx = lbox_unit * num_x
        ly = lbox_unit * num_y
        nparticles = num_x*num_y
        Jumol.create_box_by_hand(gl.molstruc,lx,ly,10.0,0.0,0.0,0.0)
        #
        num_small = convert(Int64, round(4*nparticles/(5+sqrt(5)), digits=0))
        num_large = nparticles - num_small
        type_vec = zeros(nparticles)
        type_vec[1:num_small] .= 1
        type_vec[num_small+1:length(type_vec)] .= 2
        shuffle!(type_vec)
        ## size of a unit box
        y = 0.5*lbox_unit
        atom_cntr = 1
        for i in 1:num_y
                x = 0.5*lbox_unit
                for ii in 1:num_x
                        Jumol.add_atom_by_hand(gl.molstruc,atom_cntr,type_vec[atom_cntr],x,y,0.0,0.0,0.0,0.0)
                        atom_cntr += 1
                        x += lbox_unit
                end
                y += lbox_unit
        end
end
