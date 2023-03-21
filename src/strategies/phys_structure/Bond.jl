## Bond class
######################################
######################################
######################################
# bond
#
#
#
#
#
# (c) Franz Bamer Nov-2020
######################################


## constructor of the bond
mutable struct Bond
    number::Int64
    node1
    node2
    Bond(number,node1,node2) = new(number,node1,node2)
    center_coord
end
