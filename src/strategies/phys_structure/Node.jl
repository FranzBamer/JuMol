## Node class
######################################
######################################
######################################
# node
#
#
#
#
#
# (c) Franz Bamer Nov-2020
######################################




mutable struct Node
    number::Int64
    pos_x::Float64
    pos_y::Float64
    type::Int64
    Node(number,pos_x,pos_y,type) = new(number,pos_x,pos_y,type)
end
