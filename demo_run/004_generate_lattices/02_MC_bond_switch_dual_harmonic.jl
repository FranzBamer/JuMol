## Network generation 1
######################################
######################################
######################################
# Generate a hexagonal network and
# switch manually
#
#
# (c) Franz Bamer, Nov-2020
######################################


Juno.clearconsole()


using PrettyTables
using Plots
using DelimitedFiles

## load module Jumol
include("../../src/Jumol.jl")
flush(stdout)


## input -- output
output_folder = "demo_run/004_generate_lattices/res/harmonic_dual_switch_hist/"



## run dual MC bond switch
function run_dual_MC_bond_switch(num_switches, nw)
    num_bonds = length(nw.Dualnetwork.dual_bond_list)
    cntr_sw = 0
    for i in 1:num_switches
        num_bond = rand(1:num_bonds)
        check = Jumol.switch_bond(nw,num_bond,true,true)
        if check == true
            cntr_sw += 1
            if cntr_sw % 10 == 0
                Jumol.visualize_network(nw, plot_args, true, join([output_folder,"network_",string(i),"_", string(nx),"_", string(ny),".svg"]))
                cntr_sw = 0
            end
        end
    end
end

## lattice input
nx = 15
ny = 18
network_type = "harmonic_dual"
plot_args = [false,true,false,false,false,false]
## Network object generation
nw = Jumol.Network(3.05, network_type)
Jumol.initialize_hexagonal_network(nw, nx, ny, [4,5,6,7,8,9,10], [-0.2,-0.2,-0.2,-0.2,-0.2,-0.2,-0.2])
nw.Visnetwork.dual_bond_numbering = false
nw.Visnetwork.dual_bond_textsize = 6
nw.Visnetwork.dual_node_color = "black"
nw.Visnetwork.node_color = "blue"
nw.Visnetwork.node_size = 8.0
nw.Visnetwork.node_bond_color = "blue"
nw.Visnetwork.node_bond_size = 5.0
Jumol.visualize_network(nw, plot_args, true, join([output_folder,"001_fig_hex_", string(nx),"_", string(ny),".svg"]))

## switching
run_dual_MC_bond_switch(2000, nw)
