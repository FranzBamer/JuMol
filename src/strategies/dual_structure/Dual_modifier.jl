######################################
######################################
######################################



mutable struct Dual_modifier
	dualnet
	Dual_modifier(dualnet) = new(dualnet)
end


## translate dual network
function translate_dual(mod::Dual_modifier, translate_vector)
	for node in mod.dualnet.dual_node_list
		node.pos_x += translate_vector[1]
		node.pos_y += translate_vector[2]
	end
end
