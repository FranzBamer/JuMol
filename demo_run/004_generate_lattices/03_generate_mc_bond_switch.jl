## Network generation MC
######################################
######################################
######################################
# Generate a hexagonal network and
# switch randomly with Aboav-Weaire
# save the result
# relax the silica structure
#
# (c) Franz Bamer, Dec-2020
######################################


Juno.clearconsole()


using PrettyTables
using Plots
using DelimitedFiles

## load module Jumol
include("../../src/Jumol.jl")
flush(stdout)

## calculation input variables
relax_and_save = false
# size of the network
num_x = 5
num_y = 6
bond_length= 3.05
# Aboav-Weaire parameter
alpha4 = 0.3
alpha5 = 0.3

## Network object generation and switching
nw = Jumol.Network(bond_length,"silica")

## definition of plot arguments according to the type of network
plot_args = []
if nw.network_type == "harmonic"
    plot_args = [true,true,true,true,true,false]
elseif nw.network_type == "harmonic_dual"
    plot_args = [true,true,true,true,true,false]
elseif nw.network_type == "silica"
    plot_args = [true,false,true,true,true,true]
elseif network_type == "graphene"
    plot_args = [true,false,true,true,true,true] #FIXME # not yet implemented
else
    println("cannot define plot arguments")
end

## initialize the network structure: starting with the hexagonal structure
Jumol.initialize_hexagonal_network(nw,num_x,num_y,[4,5,6,7,8,9],[0.3,0.3,0.3,0.3,0.3,0.3,0.3])
nw.Visnetwork.dual_bond_numbering = false
nw.Visnetwork.dual_bond_textsize = 6
Jumol.visualize_network(nw,plot_args)

## mc bond switch with Aboav-Weaire law (for the  first few switches no statistics)
num_switches = 10
cntr_allowed = 0
for i in 1:num_switches
    include_statistics = false
    if i > 5
        include_statistics = true
    end
    println(" ---- switch number: ",i)
    num_switch = rand(1:length(nw.Dualnetwork.dual_bond_list))
    switch_bool = Jumol.switch_bond(nw,num_switch,true,include_statistics)
    if switch_bool
        Jumol.visualize_network(nw,plot_args)
        global cntr_allowed += 1
    else
        println("this switch was not allowed")
    end
end
println("Number of switches performed: ", cntr_allowed)

## save the dual network
#Jumol.save_dual_network(nw.Dualnetwork,"output.dual")

## relax and save the output sample
if relax_and_save
    ## save the ring structure before
    Jumol.open_file(nw.Physnetwork.silicastruc.reader,"output_before_relax.lammpstrj","w")
    Jumol.write_box(nw.Physnetwork.silicastruc.reader,
                    nw.Physnetwork.silicastruc,0)
    Jumol.close_file(nw.Physnetwork.silicastruc.reader)

    ## relax the sample
    Deformer = Jumol.Aff_deform(nw.Physnetwork.silicastruc)
    num_relax = 20
    epsilon=1.0e-5
    for i in 1:num_relax
        Borelax = Jumol.Box_relax(nw.Physnetwork.silicacalc,Deformer,100,epsilon)
        Borelax.stress_tol = 1.0e-6
        Jumol.relax_biaxial(Borelax,1.0e-2,1.0e-2,1)
    end

    ## shear relax
    Borelax = Jumol.Box_relax(nw.Physnetwork.silicacalc,Deformer,100,epsilon)
    Borelax.stress_tol = 1.0e-8
    Jumol.relax_shear(Borelax, 1.0e-2, [1,2])


    ## save the ring structure after
    Jumol.open_file(nw.Physnetwork.silicastruc.reader,"output_after_relax.lammpstrj","w")
    Jumol.write_box(nw.Physnetwork.silicastruc.reader,
                    nw.Physnetwork.silicastruc,0)
    Jumol.close_file(nw.Physnetwork.silicastruc.reader)


    ## plot result
    Vis = Jumol.Vis2d(nw.Physnetwork.silicastruc)
    fig = Jumol.plot_box(Vis)
    Jumol.plot_atoms(Vis)
    display(fig)

end
