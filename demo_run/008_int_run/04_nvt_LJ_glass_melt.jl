Juno.clearconsole()
using Plots
using Random
using PrettyTables

include("../../src/Jumol.jl")

## save a function
function save_function(fname::String, vec1, vec2, factor::Int64)
        file = open(fname, "w")
        for i in 1:length(vec1)
                if i%factor == 0
                        write(file, join([string(vec1[i]), " ", string(vec2[i]), "\n"]))
                end
        end
        close(file)
end

## run the nvt script
function run_heat_bath(num_x::Int64, num_y::Int64, tau_factor::Float64)
        ## input -- output
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
        ##
        molint = Jumol.Integ(molstruc,molcalc,molunits,dt)
        tau = dt*tau_factor
        clctr = Jumol.open_file(molstruc.reader,join([folder_output,"LJ_output_nvt_LJ_melting.lammpstrj"]),"w")
        Jumol.write_box(molstruc.reader,molstruc,1)
        Jumol.run_nvt(molint, 2000, 30.0, 30.0, tau, true, 2)
        ## plot the temperature
        fig = plot(molint.timestep_vec, molint.temp_vec, color="blue", title=string(tau))
        display(fig)
        ## plot the damping parameter
        fig = plot(molint.timestep_vec, molint.zeta_vec, color="blue", title=string(tau))
        display(fig)
        ## save file
        save_function(join([folder_output,"melt_LJ_temp.dat"]), molint.timestep_vec, molint.temp_vec, 25)
        save_function(join([folder_output,"melt_LJ_zeta.dat"]), molint.timestep_vec, molint.temp_vec, 25)
end





## run  the file
#run_heat_bath(20,20,0.001)
#run_heat_bath(20,20,0.01)
#run_heat_bath(20,20,0.1)
#run_heat_bath(20,20,1.0)
#run_heat_bath(20,20,10.0)
run_heat_bath(30, 30, 100.0)
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
