## minimizer
######################################
######################################
######################################



mutable struct Units
    unit_type::Int64
    Units(unit_type=0) = new(unit_type)
    force_factor::Float64
    Kb::Float64
    barConv::Float64
end

function set_units(units::Units)

    ## unit factors
    if units.unit_type == 1
        #metal units of Lammps
        #Convert eV/A to A/ps2
        units.force_factor = 0.964853321e4

        #Boltzmann contant in Da A2/(ps2 * K)
        units.Kb = 0.831445972752 #1.38064852e-23 *0.602214076e23

        #converting bar into metal
        units.barConv = 6.022140764e-3

    else
        #Unit less
        units.force_factor = 1
        units.Kb = 1
        units.barConv = 1
    end
end
