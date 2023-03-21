Juno.clearconsole()
using Plots
using Random
using PrettyTables

include("../../src/Jumol.jl")

## run the heat bath
function run_heat_bath(num_x::Int64, num_y::Int64, tau_factor::Float64)
        ## input -- output
        folder_input = "demo_run/008_int_run/res/"
        folder_output = "demo_run/008_int_run/res/"
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
        fname = "nvt_LJ_melt_0.1.lammpstrj"
        Jumol.read_lammpstrj(molstruc, fname,1)
        #
        Jumol.initialize_structure(molstruc)
        flush(stdout)
        molcalc = Jumol.Calc(molstruc)
        Jumol.initialize_potential(molcalc,33)
        #
        ## assign masses to the particles
        for i in 1:molstruc.noa
                if molstruc.atom_list[i].type == 1
                        molstruc.atom_list[i].mass = 1.0
                elseif molstruc.atom_list[i].type == 2
                        molstruc.atom_list[i].mass = 2.0
                end
        end
        #
        ## set units
        molunits = Jumol.Units(0)
        Jumol.set_units(molunits)
        #### time integration
        dt = 1e-3
        #### melt
        molint = Jumol.Integ(molstruc,molcalc,molunits,dt)
        tau = dt*tau_factor
        fname_melt = join(["NVT_LJ_melt2_", string(tau), ".lammpstrj"])
        Jumol.run_nvt(molint, 10000, 30.0, 30.0, tau, true, 2, fname_melt)
        ## plot the temperature
        fig = plot(molint.timestep_vec, molint.temp_vec, color="blue", title=string(tau))
        display(fig)
        fname = join(["NVT_LJ_melt2_temp_", string(tau), ".pdf"])
        savefig(fig,fname)
        ## plot the damping parameter
        fig = plot(molint.timestep_vec, molint.zeta_vec, color="blue", title=string(tau))
        display(fig)
        fname = join(["NVT_LJ_melt2_zeta", string(tau), ".pdf"])
        savefig(fig,fname)
        ## save melt file
        save_function("melt2_LJ_temp.dat", molint.timestep_vec, molint.temp_vec, 25)
        save_function("melt2_LJ_zeta.dat", molint.timestep_vec, molint.zeta_vec, 25)
        #### quench the system
        molint2 = Jumol.Integ(molstruc,molcalc,molunits,dt)
        tau = dt*tau_factor
        fname_quench = join(["NVT_LJ_quench_", string(tau), ".lammpstrj"])
        Jumol.run_nvt(molint2, 2000, 30.0, 0.01, tau, true, 2, fname_quench, molint.zeta_vec[end])
        ## plot the temperature
        fig = plot(molint2.timestep_vec, molint2.temp_vec, color="blue", title=string(tau))
        display(fig)
        fname = join(["NVT_LJ_quench_temp_", string(tau), ".pdf"])
        savefig(fig,fname)
        ## plot the damping parameter
        fig = plot(molint2.timestep_vec, molint2.zeta_vec, color="blue", title=string(tau))
        display(fig)
        fname = join(["NVT_LJ_quench_zeta", string(tau), ".pdf"])
        savefig(fig,fname)
        ## save quench file
        save_function("quench_LJ_temp.dat", molint2.timestep_vec, molint2.temp_vec, 5)
        save_function("quench_LJ_zeta.dat", molint2.timestep_vec, molint2.zeta_vec, 5)
        ## plot the end result
        plot_atomic_structure(molstruc, 500, 500, 5, 8, "quenched_sample.svg",
                              [-molstruc.box.lx*0.075, molstruc.box.lx*1.25, molstruc.box.lx*1.25,-molstruc.box.lx*0.075],
                              [-molstruc.box.ly*0.075, -molstruc.box.ly*0.075, molstruc.box.ly*1.075, molstruc.box.ly*1.075])
end



## run  the file
#run_heat_bath(20,20,0.001)
#run_heat_bath(20,20,0.01)
#run_heat_bath(20,20,0.1)
#run_heat_bath(20,20,1.0)
#run_heat_bath(20,20,10.0)
run_heat_bath(20,20,100.0)
#run_heat_bath(20,20,200.0)
#run_heat_bath(20,20,300.0)
#run_heat_bath(20,20,400.0)
#run_heat_bath(20,20,500.0)
#run_heat_bath(20,20,600.0)
#run_heat_bath(20,20,700.0)
#run_heat_bath(20,20,800.0)
#run_heat_bath(20,20,900.0)
#run_heat_bath(20,20,1000.0)
#run_heat_bath(20,20,10000.0)
#run_heat_bath(20,20,100000.0)
#run_heat_bath(20,20,1000000.0)
#run_heat_bath(20,20,10000000.0)
#run_heat_bath(20,20,100000000.0)
#run_heat_bath(20,20,1000000000.0)
#run_heat_bath(20,20,10000000000.0)
#run_heat_bath(20,20,100000000000.0)
