## Atom class
######################################
######################################
######################################



import Core.Array


mutable struct Atom
	number::Int64
	type::Int64
	group::Int64
	calc_stress::Bool
	Atom(number,type,group=0,calc_stress=true) = new(number,type,group,calc_stress)
	mass::Float64
	##
	pos::Vector{Float64}
	vel::Vector{Float64}
	acc::Vector{Float64}
	##
	neighbor_indices::Vector{Int64}
	distances::Array{Float64}
	##
	pot::Float64
	force::Vector{Float64}
	## cell index, to which atom belongs to
	atom_cell_index_vec::Vector{Int64}
end
