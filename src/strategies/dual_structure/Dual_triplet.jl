## Dual triplet class
######################################
######################################
######################################
# dual triplet
#
#
#
#
#
# (c) Franz Bamer Nov-2020
######################################



mutable struct Dual_triplet
    number::Int64
    bond1
    bond2
    bond3
    Dual_triplet(number,bond1,bond2,bond3) = new(number,bond1,bond2,bond3)
    ## three nodes of the triplet
    node1
    node2
    node3
    ##
    pos_x_center::Float64
    pos_y_center::Float64
    ##
    neigh_triplet_list
end



function find_triangle_nodes(dt::Dual_triplet)
    ## the first two nodes are defined by the first bond
    dt.node1 = dt.bond1.node1
    dt.node2 = dt.bond1.node2
    ## third node from bond 2
    if dt.bond1.node1.number == dt.bond2.node1.number
        dt.node3 = dt.bond2.node2
    elseif dt.bond1.node2.number == dt.bond2.node1.number
        dt.node3 = dt.bond2.node2
    elseif dt.bond1.node2.number == dt.bond2.node2.number
        dt.node3 = dt.bond2.node1
    elseif dt.bond1.node1.number == dt.bond2.node2.number
        dt.node3 = dt.bond2.node1
    else
        println("... this should not happen ...")
    end
    ## getting the center of the triangle
    dt.pos_x_center = (dt.node1.pos_x + dt.node2.pos_x + dt.node3.pos_x)/3.0
    dt.pos_y_center = (dt.node1.pos_y + dt.node2.pos_y + dt.node3.pos_y)/3.0
end



function find_neigh_triplets(dt::Dual_triplet,triplet_list)
    ## find neighbor triplet to every bond
    dt.neigh_triplet_list = Vector()
    ## bond 1
    bond_num_list = [dt.bond1.number dt.bond2.number dt.bond3.number]
    for i in 1:3 # every triplet has three bonds
        bond_num = bond_num_list[i]
        cntr = 1
        for ii in 1:length(triplet_list)
            triplet_candidate = triplet_list[ii]
            if dt.number != triplet_candidate.number
                if bond_num == triplet_candidate.bond1.number || bond_num == triplet_candidate.bond2.number || bond_num == triplet_candidate.bond3.number
                    push!(dt.neigh_triplet_list,cntr)
                    break
                end
            end
            cntr += 1
        end
    end
    #println("neighbor triplet list of triplet ", dt.number)
    #println(dt.neigh_triplet_list)
end

## find adjacent ring triplets -> this function is used to proceed from triplet to triplet around a node
## input: Dual_triplet (self) , ring_triplet_list: is the list of triplets defining a ring
function find_adjacent_ring_triplets(dt::Dual_triplet,ring_triplet_list)
    adjacent_ring_triplets = Vector()
    for i in 1:length(dt.neigh_triplet_list)
        neigh_triplet_num = dt.neigh_triplet_list[i]
        for ii in 1:length(ring_triplet_list)
            ring_triplet_num = ring_triplet_list[ii]
            if neigh_triplet_num == ring_triplet_num
                push!(adjacent_ring_triplets,ring_triplet_num)
            end
        end
    end
    return adjacent_ring_triplets
end

## find third node of the triplet
function find_third_node_of_triplet(dt::Dual_triplet,bond)
    tr_node_list = [dt.node1 dt.node2 dt.node3]
    bond = identify_bond_by_number(dt, bond.number)
    third_num = 0
    for node_num in tr_node_list
        if node_num != bond.node1 && node_num != bond.node2
            third_num = node_num
        end
    end
    return third_num
end

## internal function -> identify a bond by its number
function identify_bond_by_number(dt::Dual_triplet,bond_number)
    if bond_number == dt.bond1.number
        return dt.bond1
    end
    if bond_number == dt.bond2.number
        return dt.bond2
    end
    if bond_number == dt.bond3.number
        return dt.bond3
    end
end
