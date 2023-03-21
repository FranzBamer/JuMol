######################################
######################################
######################################
## set the box and the atoms to affine deformations
# four atom benchmark example

Juno.clearconsole()


using PrettyTables
using Plots


include("../../src/Jumol.jl")
flush(stdout)


molstruc = Jumol.Structure()
molstruc.rc = 2.0
molstruc.rskin = 0.4
molstruc.pbx = 1
molstruc.pby = 1

Jumol.initialize_structure_objects(molstruc)
Jumol.create_box_by_hand(molstruc,5.0,5.0,10.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,1,1,1.0,1.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,2,1,3.0,2.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,3,1,2.0,4.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,4,1,4.0,4.0,0.0,0.0,0.0,0.0)
Jumol.initialize_structure(molstruc)

## initial box dimeions
lx0 = molstruc.box.lx
ly0 = molstruc.box.ly

Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = 1000
Plotter.hfig = 500
fig = Jumol.plot_box(Plotter)
Jumol.plot_atoms(Plotter,"blue",5)
display(fig)

## plot initial structure
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-lx0*0.075, lx0*2.0, lx0*2.0,-lx0*0.075],
                      [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)

## load deformation object
Deformer = Jumol.Aff_deform(molstruc)
Jumol.set_affine_deform_shear(Deformer,4.0,0.0,0.0)

## plot atomic structure deformed
Plotter = Jumol.Vis2d(molstruc)
Plotter.bfig = bfig
Plotter.hfig = hfig
fig = Jumol.plot_box(Plotter)
Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
                      [-lx0*0.075, lx0*2.0, lx0*2.0,-lx0*0.075],
                      [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
display(fig)

## done
println("...Done...")
