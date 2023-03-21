## run a minimization of one configuration
######################################
######################################
######################################

Juno.clearconsole()


#using Revise
using PrettyTables
using Plots
using Random

## load module Jumol
include("../../src/Jumol.jl")
flush(stdout)


## plot functions
function plot_atomic_structure(molstruc, hfig::Int64, bfig::Int64, size1::Int64, size2::Int64, filename::String,
                               framex::Vector, framey::Vector)
    Plotter = Jumol.Vis2d(molstruc)
    Plotter.bfig = bfig
    Plotter.hfig = hfig
    fig = Jumol.plot_box(Plotter)
    Jumol.plot_white_frame(Plotter,framex, framey)
    Jumol.plot_atoms_type(Plotter, 1, "blue", "lightgray", size1)
    Jumol.plot_atoms_type(Plotter, 2, "blue", "lightblue", size2)
    savefig(filename)
    display(fig)
end




function generate_LJ_glass(num::Int64, num_x::Int64, num_y::Int64, size1::Float64, size2::Float64, relax::Bool, savefile::Bool)
        ## create Jumol object
        molstruc = Jumol.Structure()
        molstruc.rc = 2.4
        molstruc.rskin = 0.4
        molstruc.pbx = 1
        molstruc.pby = 1
        #
        Jumol.initialize_structure_objects(molstruc)
        #
        ## create lj lattice
        Gen_Lat = Jumol.Gen_lattice(molstruc)
        Jumol.create_random_pos_LJ(Gen_Lat, num_x, num_y)
        #
        Jumol.initialize_structure(molstruc)
        flush(stdout)
        molcalc = Jumol.Calc(molstruc)
        Jumol.initialize_potential(molcalc,33)
        #
        ## assign masses to the particles
        for i in 1:molstruc.noa
                molstruc.atom_list[i].mass = 1.0
        end
        #
        ## set units
        molunits = Jumol.Units(0)
        Jumol.set_units(molunits)
        ## plot atomic structure initial
        Plotter = Jumol.Vis2d(molstruc)
        Plotter.bfig = bfig
        Plotter.hfig = hfig
        fig = Jumol.plot_box(Plotter)
        Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                              [-molstruc.box.lx*0.075, molstruc.box.lx*1.075, molstruc.box.lx*1.075,-molstruc.box.lx*0.075],
                              [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075] )
        display(fig)
        #
        ## run minimization
        Minimizer = Jumol.Min(molstruc,molcalc)
        Minimizer.alpha_min = 1.0e0
        Jumol.run_cg(Minimizer, 100000, 1.0e-6)
        Jumol.put_atoms_back_to_box(molstruc)
        #
        if relax == true
                ## deformation object
                Deformer = Jumol.Aff_deform(molstruc)
                ## athermal box relaxation
                ## athermal shear relaxation
                Borelax = Jumol.Box_relax(molcalc,Deformer,100,1.0e-6)
                Jumol.relax_vol(Borelax, 1.0e-1)
                #Jumol.relax_biaxial(Borelax, 1.0e-1, 1.0e-1)
                Jumol.relax_shear(Borelax, 5.0e-2, [1,2])
        end
        #
        ## plot atomic structure again
        Plotter = Jumol.Vis2d(molstruc)
        Plotter.bfig = bfig
        Plotter.hfig = hfig
        fig = Jumol.plot_box(Plotter)
        Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                              [-molstruc.box.lx*0.075, molstruc.box.lx*1.075, molstruc.box.lx*1.075,-molstruc.box.lx*0.075],
                              [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075] )
        display(fig)
        #
        if savefile == true
                output_folder = "demo_run/004_generate_lattices/res/LJ_glass/"
                fname = join([output_folder,"sample_",string(num_x),"_",string(num_y),"_num_",string(num),".lammpstrj"])
                Jumol.open_file(molstruc.reader, fname, "w")
                Jumol.write_box(molstruc.reader, molstruc, 1)
                Jumol.close_file(molstruc.reader)
        end
        #
        println("...LJ glass generated...")
end

#### ball sizes: small, large
# samples 10 x 10 --> 12.5, 22.5
x_num = 10
y_num = 10
num_samples = 1
sample_start = 1
size1 = 12.5
size2 = 22.5
for i in sample_start:sample_start
        generate_LJ_glass(i, x_num, y_num, size1, size2, true, true)
end
