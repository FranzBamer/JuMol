
Juno.clearconsole()
using Plots

include("../../src/Jumol.jl")

## plot properties
size1 = 13.5
size2 = 25.0
hfig = 500
bfig = 680

## structure
molstruc = Jumol.Structure()
molstruc.rc = 2.4
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
## assign masses to the particles
for i in 1:molstruc.noa
        molstruc.atom_list[i].mass = 10.0
end
flush(stdout)

## input -- output
folder_output = "demo_run/008_int_run/res/"

## calculation object
molcalc = Jumol.Calc(molstruc)
## including Leonard Jones potential
Jumol.initialize_potential(molcalc,33)
molunits = Jumol.Units(1)
Jumol.set_units(molunits)
#
clctr = Jumol.open_file(molstruc.reader,join([folder_output,"LJ_output_nvt_4_atoms.lammpstrj"]),"w")
Jumol.write_box(molstruc.reader,molstruc,1)
#
molint = Jumol.Integ(molstruc,molcalc,molunits,1.0e-3)
Jumol.run_nvt(molint,1000,1000,10000,2.0)

## Done
println("..Done..")
