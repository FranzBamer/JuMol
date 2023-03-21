
Juno.clearconsole()
using Plots

include("../../src/Jumol.jl")

## input -- output
folder_input = "demo_run/999_benchmark_samples/2D_silica_8x8/sig_1.00/"
folder_output = "demo_run/008_int_run/res/"

molstruc = Jumol.Structure()
molstruc.rc = 2.4
molstruc.rskin = 0.4
molstruc.pbx = 1
molstruc.pby = 1
#
Jumol.initialize_structure_objects(molstruc)
Jumol.read_lammpstrj(molstruc,join([folder_input,"sample_0.lammpstrj"]))
#
Jumol.initialize_structure(molstruc)
#
for i in 1:molstruc.noa
    if molstruc.atom_list[i].type == 1
        molstruc.atom_list[i].mass = 1.0
    else
        molstruc.atom_list[i].mass = 15.999/28.0855
    end
end
flush(stdout)
molcalc = Jumol.Calc(molstruc)
Jumol.initialize_potential(molcalc,13)
Jumol.calc_all_pair_forces(molcalc)
molunits = Jumol.Units(0)
Jumol.set_units(molunits)
#
clctr = Jumol.open_file(molstruc.reader,join([folder_output,"network_glass_nvt.lammpstrj"]),"w")
Jumol.write_box(molstruc.reader,molstruc,1)
molint = Jumol.Integ(molstruc,molcalc,molunits,1.0e-3)
Jumol.run_nvt(molint,10,0.01,10.0,1.0)
Jumol.close_file(molstruc.reader)

## Done
println("..Done..")
