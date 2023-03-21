## Box class
######################################
######################################
######################################


mutable struct Box
	Box() = new()
	## box dimensions
	lx::Float64
	ly::Float64
	lz::Float64
	lxy::Float64
	lyz::Float64
	lxz::Float64
	## basis vectors
	h1::Vector{Float64}
	h2::Vector{Float64}
	h3::Vector{Float64}
	## length of the tricilic box
	l1::Float64
	l2::Float64
	l3::Float64
	## unit vectors
	e1::Vector{Float64}
	e2::Vector{Float64}
	e3::Vector{Float64}
end



function set_box_basis_vectors(box::Box)
	## box basis vectors
	box.h1 = [box.lx,0,0]
	box.h2 = [box.lxy,box.ly,0]
	box.h3 = [box.lxz,box.lyz,box.lz]
	## lenght of the triclinic box
	box.l1 = norm(box.h1)
	box.l2 = norm(box.h2)
	box.l3 = norm(box.h3)
	## box unit vectors
	box.e1 = box.h1/box.l1
	box.e2 = box.h2/box.l2
	box.e3 = box.h3/box.l3
end
