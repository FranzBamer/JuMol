## Cell class
######################################
######################################
######################################

import Core.Array


mutable struct Cell

	## cell dimensions
	cell_id::Int64
	Cell(cell_id, atoms_in_cell_list=Vector()) = new(cell_id, atoms_in_cell_list)

	atoms_in_cell_list::Vector{Int64}

end
