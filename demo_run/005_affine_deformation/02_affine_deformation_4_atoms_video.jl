######################################
######################################
######################################
## affine displacements video
# four atom benchmark example

Juno.clearconsole()


include("../../src/Jumol.jl")
flush(stdout)


molstruc = Jumol.Structure()

molstruc.rc = 2.0
molstruc.rskin = 0.4
molstruc.pbx = 1
molstruc.pby = 1

## input -- output
output_folder = "demo_run/005_affine_deformation/res/"

Jumol.initialize_structure_objects(molstruc)
Jumol.create_box_by_hand(molstruc,5.0,5.0,10.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,1,1,1.0,1.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,2,1,3.0,2.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,3,1,2.0,4.0,0.0,0.0,0.0,0.0)
Jumol.add_atom_by_hand(molstruc,4,1,4.0,4.0,0.0,0.0,0.0,0.0)
Jumol.initialize_structure(molstruc)

## initial box_dimensions
lx0 = molstruc.box.lx
ly0 = molstruc.box.ly

## Deformation object
Deformer = Jumol.Aff_deform(molstruc)

## run affine deformation history
n = 100
delta_lxy = 0.01
anim = @animate for i in 1:n
	println("Deformation step: ",i)
	Jumol.set_affine_deform_shear(Deformer,delta_lxy,0.0,0.0)
	## plot atomic structure deformed
	Plotter = Jumol.Vis2d(molstruc)
	Plotter.bfig = bfig
	Plotter.hfig = hfig
	fig = Jumol.plot_box(Plotter)
	Jumol.plot_atomic_structure_binary(Plotter, size1, size2,
	                      [-lx0*0.075, lx0*1.50, lx0*1.50,-lx0*0.075],
	                      [-ly0*0.075, -ly0*0.075, ly0*1.075, ly0*1.075] )
	display(fig)
end
gif(anim, join([output_folder,"benchmark_video_4atoms.gif"]), fps = 25)



println("...Done...")
