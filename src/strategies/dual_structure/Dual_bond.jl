## Dual bond class
######################################
######################################
######################################
# dual bond
#
#
#
#
#
# (c) Franz Bamer Nov-2020
######################################



mutable struct Dual_bond
    number::Int64
    node1
    node2
    Dual_bond(number,node1,node2) = new(number,node1,node2)
    ##
    neigh_triplets
end

## find the two  neighboring triangles
function find_neighboring_triplets(bd::Dual_bond,triplet_list)
    bd.neigh_triplets = Vector()
    for i in 1:length(triplet_list)
        tr = triplet_list[i]
        if tr.bond1.number == bd.number || tr.bond2.number == bd.number || tr.bond3.number == bd.number
            push!(bd.neigh_triplets,tr)
        end
    end
    return bd.neigh_triplets
end
