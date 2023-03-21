
Juno.clearconsole()

using Plots
include("../../src/Jumol.jl")

## plot properties
size1 = 13.5
size2 = 25.0
hfig = 500
bfig = 680

## input -- output
folder_output = "demo_run/008_int_run/res/"

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
flush(stdout)



## assign masses to the particles
for i in 1:molstruc.noa
        molstruc.atom_list[i].mass = 28.0855
end
molcalc = Jumol.Calc(molstruc)

##including Potential
Jumol.initialize_potential(molcalc,33)
molunits = Jumol.Units(1)
Jumol.set_units(molunits)
clctr = Jumol.open_file(molstruc.reader,join([folder_output,"output_nve.lammpstrj"]),"w")
Jumol.write_box(molstruc.reader,molstruc,1)
molint = Jumol.Integ(molstruc,molcalc,molunits,0.001)
Jumol.run_nve(molint,1000)
Jumol.write_box(molstruc.reader,molstruc,2)
Jumol.close_file(molstruc.reader)


println("..Done..")
