## Dual node class
######################################
######################################
######################################
# dual node
#
#
#
#
#
# (c) Franz Bamer Nov-2020
######################################



mutable struct Dual_node
    number::Int64
    type::Int64
    pos_x::Float64
    pos_y::Float64
    Dual_node(number,type,pos_x,pos_y) = new(number,type,pos_x,pos_y)
    adjacent_node_list
    adjacent_triplet_list
	adjacent_triplet_list_ordered #ring
    ## for Aboav-Weaire law
    mean_adjacent_node_type::Float64 #(mean adjacent ring size of the physical network)
    alpha_aboav_weaire::Float64
end

## get the adjacent neighbor nodes
function get_adjacent_neigh_nodes(nd::Dual_node,dual_node_list,dual_bond_list)
	nd.adjacent_node_list = Vector()
	for i in 1:length(dual_bond_list)
		bond = dual_bond_list[i]
		node_num_bond_list = [bond.node1.number bond.node2.number]
		## is the node part of this bond?
		bond_check = false
		for ii in 1:length(node_num_bond_list)
			node_num = node_num_bond_list[ii]
			if node_num == nd.number
				bond_check = true
			end
		end
		## find the adjacent node if the node is part of this bond
		if bond_check==true
			for ii in 1:length(node_num_bond_list)
				node_num = node_num_bond_list[ii]
				if node_num != nd.number
					push!(nd.adjacent_node_list,node_num)
				end
			end
		end
	end
	## calculate the mean surrounding node type (for Aboav-Weaire law)
	nd.mean_adjacent_node_type = 0.0
	for i in 1:length(nd.adjacent_node_list)
		node_num = nd.adjacent_node_list[i]
		nd.mean_adjacent_node_type += float(dual_node_list[node_num].type)
	end
	nd.mean_adjacent_node_type = nd.mean_adjacent_node_type/float(nd.type)
end

## find neighboring triplets to the node in an ordered manner (ring)
## every triplet must know his three neighbors
function find_adjacent_neigh_triplets_ordered(nd::Dual_node,dual_triplet_list)
	nd.adjacent_triplet_list = Vector()
	for i in 1:length(dual_triplet_list)
		tr = dual_triplet_list[i]
		n1 = tr.node1.number
		n2 = tr.node2.number
		n3 = tr.node3.number
		node_list = [n1 n2 n3]
		for ii in 1:length(node_list)
			node_num = node_list[ii]
			if node_num == nd.number
				push!(nd.adjacent_triplet_list,tr.number)
			end
		end
	end
	## order the triplet list
	nd.adjacent_triplet_list_ordered = Vector()
	## define first triplet
    triplet1_num = nd.adjacent_triplet_list[1]
    ## find the three neighbor triplets around the triplet
    find_neigh_triplets(dual_triplet_list[triplet1_num],dual_triplet_list)
    ## find those two neighbor triplets that belong to the ring
    adjacent_ring_triplets = find_adjacent_ring_triplets(dual_triplet_list[triplet1_num],nd.adjacent_triplet_list)
    ## define the starting direction by the first two adjacent triplets randomly
    triplet2_num = adjacent_ring_triplets[1]
    push!(nd.adjacent_triplet_list_ordered,triplet1_num)
    push!(nd.adjacent_triplet_list_ordered,triplet2_num)
	for i in 2:length(nd.adjacent_triplet_list)-1
        ## find the three neighbor triplets around the next triplet
        triplet_num = nd.adjacent_triplet_list_ordered[i]
        find_neigh_triplets(dual_triplet_list[triplet_num],dual_triplet_list)
        ## find those two neighbor triplets that belong to the ring
        adjacent_ring_triplets = find_adjacent_ring_triplets(dual_triplet_list[triplet_num],nd.adjacent_triplet_list)
		## go forward
        next_triplet_num = 0
		for ii in 1:length(adjacent_ring_triplets)
			next_num = adjacent_ring_triplets[ii]
			num_before = nd.adjacent_triplet_list_ordered[i-1]
			if next_num != num_before
				next_triplet_num = next_num
			end
		end
        ## add triplet to ring_list
        push!(nd.adjacent_triplet_list_ordered, next_triplet_num)
	end
	#println("ring associated to dual node: ", nd.number)
	#println(nd.adjacent_triplet_list_ordered)
end
