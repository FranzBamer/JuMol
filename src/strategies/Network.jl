## Basic network class
######################################
######################################
######################################
# dual network
# physical network
# Visualization class
# Ring statistics class
#
# Network types:   - harmonic
#                  - harmonic dual
#                  - silica yukawa
#                  - graphene
#
# (c) Franz Bamer March-2022
######################################

include("dual_structure/Dual_network.jl")
include("phys_structure/Phys_network.jl")
include("dual_structure/Ring_statistics.jl")
include("../vis/Network_vis.jl")

## constructor of the Network class
mutable struct Network
    bond_length::Float64
	network_type::String
	pbc
    Network(bond_length, network_type, pbc=true) = new(bond_length, network_type, pbc)
    Dualnetwork
    Physnetwork
	Visnetwork
	ring_size_list_allowed
	ring_stat_vec
	## ring statistics object
	Ring_stat
	##
	cntr_attempted_switch::Float64
	cntr_top_accepted_switch::Int64
	cntr_switch::Int64
	##
	full_min_rate_factor::Int64
end

## initialize hexagonal network
function initialize_hexagonal_network(nw::Network, num_x, num_y, sig,
	                                  ring_size_list_allowed=[4,5,6,7,8,9,10],
									  alpha_target_list = [0.3,0.3,0.3,0.3,0.3,0.3,0.3],
									  target_stat = [0.0399705727098,0.275646469643,0.401653587686,0.212793419493,0.0578914101397,0.0101351126893,0.00132376850128],
									  full_min_rate_factor=1,
									  include_lammps=true, lammps_executable="None")
	##
	nw.full_min_rate_factor = full_min_rate_factor
	## define allowed ring sizes of the network
	nw.ring_size_list_allowed = ring_size_list_allowed
	## define ring statistics
	nw.ring_stat_vec = target_stat
	println(" -- initializing hexagonal network -- ")
	## check if num_y is even: if not stop the script!
	if num_y%2 != 0
		throw(DomainError(num_y, "The replication of dual rows (num_y) must be an even number"))
	end
    ## initialize the dual network
    nw.Dualnetwork = Dual_network(nw.bond_length)
    build_starting_dual_hex(nw.Dualnetwork,num_x,num_y,nw.ring_size_list_allowed)
    ## create the physical equivalent
    nw.Physnetwork = Phys_network(nw.Dualnetwork)
	nw.Physnetwork.include_lammps = true
	nw.Physnetwork.lammps_executable = lammps_executable
	nw.Physnetwork.sig = sig
    # load the dual into the physical network
    dual2phys(nw.Physnetwork)
	if nw.network_type == "harmonic"
		build_phys_molstruc(nw.Physnetwork)
	elseif nw.network_type == "harmonic_dual"
		build_phys_molstruc(nw.Physnetwork)
	elseif nw.network_type == "silica"
		build_phys_molstruc(nw.Physnetwork)
		build_silica_molstruc(nw.Physnetwork)
	else
		println("Type of network not recognized!")
	end
	## initialize the visualization object
	nw.Visnetwork = Network_vis(nw.Dualnetwork,nw.Physnetwork)
	initialize_standard_plot_settings(nw.Visnetwork)
	## initialize ring statistics object
	nw.Ring_stat = Ring_statistics(nw.Dualnetwork, nw.ring_size_list_allowed, target_stat)
	calculate_all(nw.Ring_stat)
	assign_prev_stat(nw.Ring_stat)
	nw.Ring_stat.alpha_aboave_weaire_list_target = alpha_target_list
	plot_all(nw.Ring_stat)
	## set all counters
	nw.cntr_attempted_switch = 0
	nw.cntr_top_accepted_switch = 0
	nw.cntr_switch = 0
end

## switch a bond
function switch_bond(nw::Network,bond_number,plot_info=false,consider_ring_stat=false)
	##
	nw.cntr_attempted_switch += 1
	## general boolean of switch
	switch_bool = true
	## check the objective function of the ring statistics before the switch
	chi_before = calc_chi(nw.Ring_stat)
	## switch the dual bond
	println("switching bond: ",bond_number)
	check_ring_sizes = switch_bond(nw.Dualnetwork,bond_number)
	if check_ring_sizes
		## calculating ring statistics
		calculate_all(nw.Ring_stat)
		chi_after = calc_chi(nw.Ring_stat)
		## step downwards or upwards?
		delta_chi = chi_after-chi_before
		println("delta_chi: "*string(delta_chi))
		check_stat = true
		if consider_ring_stat
			check_stat = make_decision(nw.Ring_stat, delta_chi, 1.0e-4)
		end
		if plot_info
			plot_all(nw.Ring_stat)
		end
		if check_stat # statistical change allowed (last level, switch finally accepted)
			##
			nw.cntr_top_accepted_switch += 1
			## create the physical equivalent
			println("- switch allowed according to Aboav-Weaire -")
			coordination_check = true
			## in case of a silica network
			if nw.network_type == "harmonic"
				assign_prev_stat(nw.Ring_stat)
				dual2phys(nw.Physnetwork)
				build_phys_molstruc(nw.Physnetwork)
				minimize_phys_network(nw.Physnetwork)
				update_node_pos(nw.Physnetwork)
				phys2dual(nw.Dualnetwork,nw.Physnetwork)
				switch_bool = true
				nw.cntr_switch += 1
				chi_before = chi_after
			elseif nw.network_type == "harmonic_dual"
				assign_prev_stat(nw.Ring_stat)
				build_dual_molstruc(nw.Dualnetwork)
				minimize_dual_molstruc(nw.Dualnetwork)
				update_dual_network(nw.Dualnetwork)
				dual2phys(nw.Physnetwork)
				build_phys_molstruc(nw.Physnetwork)
				nw.cntr_switch += 1
				swith_bool = true
				chi_before = chi_after
			elseif nw.network_type == "silica"
				dual2phys(nw.Physnetwork)
				build_phys_molstruc(nw.Physnetwork)
				minimize_phys_network(nw.Physnetwork)
				update_node_pos(nw.Physnetwork)
				if nw.cntr_top_accepted_switch % nw.full_min_rate_factor == 0
					build_silica_molstruc(nw.Physnetwork)
					minimize_silica_network(nw.Physnetwork)
					## silica consistency check
					coordination_check = silica_consistency_check(nw.Physnetwork,2,1.8)
					if coordination_check
						assign_prev_stat(nw.Ring_stat)
						update_silica_node_pos(nw.Physnetwork)
						phys2dual(nw.Dualnetwork,nw.Physnetwork)
						nw.cntr_switch += nw.full_min_rate_factor
						chi_before = chi_after
					else
						## switch back
						switch_bond(nw.Dualnetwork,bond_number)
						println("switch not allowed, coordination check not fullfilled!")
						println("switching back...")
						switch_bool = false
					end
				end
			elseif nw.network_type == "graphene"
				coordination_check = false #FIXME # graphene must be implemented
				if coordinateion_check
					# FIXME
				else
					# FIXME
				end
			else
				println("network type not recognized!")
			end
		else # statistical change not allowed
			## switch back
			println("switch not allowed, ring statistics not allowed!")
			println("switching back...")
			switch_bond(nw.Dualnetwork,bond_number)
			# recalculating ring statistics
			calculate_all(nw.Ring_stat)
			switch_bool = false
		end
	else
		## switch back
		println("switch not allowed, one or more ring sizes are not in the ring list!")
		println("switching back...")
		switch_bond(nw.Dualnetwork,bond_number)
		switch_bool = false
	end
	return switch_bool
end

## switch a set of bonds on the dual lattice and minimize once
#  using the harmonic potential and then the silica potential
function switch_bond_set(nw::Network, set_of_bond_numbers,
	                     plot_info=false, consider_ring_stat=false,
						 include_lammps=true)
	##
	nw.cntr_attempted_switch += length(set_of_bond_numbers)
	## general boolean of switch
	switch_bool = true
	## check the objective function of the ring statistics before the switch
	calculate_all(nw.Ring_stat)
	chi_before = calc_chi(nw.Ring_stat)
	## create a list of the actually performed dual switches
	switch_list = Vector()
	chi_list_alpha = Vector()
	chi_list_het = Vector()
	## run the switch attempts
	for i in 1:length(set_of_bond_numbers)
		## bond number to switch
		bond_number = set_of_bond_numbers[i]
		## switch the dual bond
		println("switching bond: ",bond_number)
		## first check: are the sizes of rings involved allowed?
		check_ring_sizes = switch_bond(nw.Dualnetwork, bond_number)
		if check_ring_sizes
			check_stat = true
			if consider_ring_stat
				## second check: is it downhill according to the objective function chi = chi_het + chi_alpha
				calculate_all(nw.Ring_stat)
				chi_after = calc_chi(nw.Ring_stat)
				delta_chi = chi_after - chi_before
				println("delta_chi: "*string(delta_chi))
				check_stat = make_decision(nw.Ring_stat,delta_chi,1.0e-4)
				if check_stat
					println("Metropolis decision: accept")
				else
					println("Metropolis decision: reject")
				end
				#plot_all(nw.Ring_stat)
				#readline()
			end
			if check_stat
				## switch is accepted according to the ring statistics
				push!(switch_list, bond_number)
				push!(chi_list_alpha, nw.Ring_stat.chi_alpha)
				push!(chi_list_het, nw.Ring_stat.chi_het)
				## assign new chi value for the next potential switch
				chi_before = chi_after
			else
				## switch is not accepted according to the ring statistics
				switch_bond(nw.Dualnetwork, bond_number)
			end
		else #check_ring_sizes == false
			# switch back
			switch_bond(nw.Dualnetwork, bond_number)
		end
	end
	println("accepted switch_list...")
	println(switch_list)
	#readline()
	if length(switch_list) > 0
		## minimize harmonic potential
		dual2phys(nw.Physnetwork)
		build_phys_molstruc(nw.Physnetwork)
		minimize_phys_network(nw.Physnetwork)
		update_node_pos(nw.Physnetwork)
		## minimize silica potential
		coordination_check = false
		build_silica_molstruc(nw.Physnetwork)
		if include_lammps
			coordination_check = minimize_silica_network_lammps_plus_coord_check(nw.Physnetwork)
		else
			minimize_silica_network(nw.Physnetwork, 1.001)
			coordination_check = silica_consistency_check(nw.Physnetwork,2,1.8)
		end
		if coordination_check
			update_silica_node_pos(nw.Physnetwork)
			phys2dual(nw.Dualnetwork,nw.Physnetwork)
			nw.cntr_switch += length(switch_list)
		else
			## swap switch_list
			reverse_switch_list = swap_list(nw, switch_list)
			#println("reverse_switch_list")
			#println(reverse_switch_list)
			#readline()
			## switch back
			println("switch not allowed, coordination check not fullfilled!")
			println("switching back...")
			for num_bond in reverse_switch_list
				switch_bond(nw.Dualnetwork,num_bond)
			end
			switch_bool = false
		end
	else
		switch_bool = false
	end
	return switch_bool, switch_list, chi_list_alpha, chi_list_het
end

## visualize the network structure
# arg 1: network object
# arg 2: plot arguments in a list of booleans
# entry 1: plot nodes of the physical (harmonic) network (boolean)
# entry 2: plot bonds of the physical (harmonic) network (boolean)
# entry 3: plot triplets of the dual netork (boolean)
# entry 4: plot bonds of the dual (boolean)
# entry 5: plot nodes of the dual (boolean)
# entry 6: plot the silica network (boolean)
function visualize_network(nw::Network, plot_args=[true,true,false,false,false,false], save="false", figname="fig_output.svg")
	fig = plot_box(nw.Visnetwork)
	if plot_args[1] #physical node plots
		plot_physical_nodes(nw.Visnetwork)
	end
	if plot_args[2] #physical_bond_plot
		plot_physical_bonds(nw.Visnetwork)
	end
	if plot_args[3] #dual_triplet_plot
		plot_dual_triplets(nw.Visnetwork)
	end
	if plot_args[4] #dual_bond_plot
		plot_bonds(nw.Visnetwork)
	end
	if plot_args[5] #dual_node_plot
		plot_dual_nodes(nw.Visnetwork)
	end
	if plot_args[6] #silica_network_plot
		plot_physical_silica_bonds(nw.Visnetwork)
	end
	display(fig)
	if save == true
		savefig(fig,figname)
	end
	return fig
end



## swap a list
function swap_list(nw::Network, input_list)
	output_list = Vector()
	len = length(input_list)
	for i in 1:len
		push!(output_list,input_list[len+1-i])
	end
	return output_list
end
