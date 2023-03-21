######################################
######################################
######################################
## affine displacements video
# four atom benchmark example

Juno.clearconsole()


include("../../src/Jumol.jl")
flush(stdout)


molstruc = Jumol.Structure()
molstruc.rc = 10.0
molstruc.rskin = 0.5
molstruc.pbx = 1
molstruc.pby = 1

## input -- output
input_folder = "demo_run/999_benchmark_samples/2D_silica_8x8/sig_1.00/"
output_folder = "demo_run/005_affine_deformation/res/"

Jumol.initialize_structure_objects(molstruc)
filename = join([input_folder,"sample_0.lammpstrj"])
Jumol.read_lammpstrj(molstruc,filename)
Jumol.initialize_structure(molstruc)

## initial box_dimensions
lx0 = molstruc.box.lx
ly0 = molstruc.box.ly

## plot properties
size1 = 2.75
size2 = 3.5
hfig = 500
bfig = 680

## run affine deformation
Deformer = Jumol.Aff_deform(molstruc)
n = 100
delta_lxy = 0.25
anim = @animate for i in 1:n
	println("deformation step: ",i)
	Jumol.set_affine_deform_shear(Deformer,delta_lxy,0.0,0.0)
	#Jumol.set_affine_deform_vol(molstruc,5.0,0.0,0.0)
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
gif(anim, join([output_folder,"benchmark_video_2D_silica.gif"]), fps = 25)



println("...Done...")
