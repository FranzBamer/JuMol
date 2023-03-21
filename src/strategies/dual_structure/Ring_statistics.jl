## Ring statistics class
######################################
######################################
######################################
# ring histogram
# mean and variance
# pi_matrix
# Aboav-Weaire parameter
#
#
# (c) Franz Bamer Nov-2020
######################################

using PrettyTables

## constructor of ring statistics
mutable struct Ring_statistics
    Dualnetwork
    r_list
	target_stat_vec
    Ring_statistics(Dualnetwork,r_list,target_stat_vec) = new(Dualnetwork,r_list,target_stat_vec)
	##
	ring_count_mat::Array{Int64}
	ring_stat_mat::Array{Float64}
    mu::Float64
    sig::Float64
    pi_matrix::Array{Float64}
	alpha_aboav_weaire_list::Vector{Float64}
	##
	ring_count_mat_prev::Array{Int64}
	ring_stat_mat_prev::Array{Float64}
    mu_prev::Float64
    sig_prev::Float64
    pi_matrix_prev::Array{Float64}
	alpha_aboav_weaire_list_prev::Vector{Float64}
	## target parameters
	alpha_aboave_weaire_list_target::Vector{Float64}
	## objective function
	chi_het::Float64
	chi_alpha::Float64
end

## assign previous ring statistics for future switches
function assign_prev_stat(rs::Ring_statistics)
	rs.ring_count_mat_prev = copy(rs.ring_count_mat)
	rs.ring_stat_mat_prev = copy(rs.ring_stat_mat)
	rs.mu_prev = copy(rs.mu)
	rs.sig_prev = copy(rs.sig)
	rs.pi_matrix_prev = copy(rs.pi_matrix)
	rs.alpha_aboav_weaire_list_prev = copy(rs.alpha_aboav_weaire_list)
end

## calc all statistics
function calculate_all(rs::Ring_statistics)
	calc_hist(rs)
	calc_pi_matrix(rs)
	calc_aboav_weaire(rs)
end

## plot all statistical information
function plot_all(rs::Ring_statistics)
	println("#######################################")
	println("########### Ring statistics ###########")
	println("Allowed ring sizes:")
	println(rs.r_list)
	ring_hist_compare_mat = zeros(size(rs.ring_count_mat,1)+2,size(rs.ring_count_mat,2))
	ring_hist_compare_mat[1:size(rs.ring_count_mat,1),1:size(rs.ring_count_mat,2)] = rs.ring_count_mat
	ring_hist_compare_mat[3,:] = rs.ring_stat_mat[2,:]
	ring_hist_compare_mat[4,:] = rs.target_stat_vec
	println("Ring count histogram:")
	pretty_table(ring_hist_compare_mat,
		         noheader = true, crop = :horizontal)
	println("Ring stat histogram:")
	pretty_table(rs.ring_stat_mat,
			     noheader = true, crop = :horizontal,
				 formatters = ft_round(3))
	println("Mean ring size: ", rs.mu)
	println("Standard deviation ring size: ", rs.sig)
	println("Aboav-Weaire parameters: ")
	pretty_table(rs.alpha_aboav_weaire_list',
			     noheader = true, crop = :horizontal,
				 formatters = ft_round(3))
	println("PI-matrix:")
	pretty_table(rs.pi_matrix,
		         noheader = true, crop = :horizontal,
		         formatters = ft_round(3))
end

## calc histogram
function calc_hist(rs::Ring_statistics)
	## ring statistics matrix
	rs.ring_count_mat = zeros(2,20)
	# fill first column (n-fold rings)
	for i = 1:size(rs.ring_count_mat,2)
		rs.ring_count_mat[1,i] = i
	end
	# fill second column (ring counts)
	for dn in rs.Dualnetwork.dual_node_list
		rs.ring_count_mat[2,dn.type] += 1
	end
	# create ring stat matrix (area of the histogram is one)
	rs.ring_stat_mat = zeros(2,20)
	for i in 1:size(rs.ring_stat_mat,2)
		rs.ring_stat_mat[1,i] = i
		rs.ring_stat_mat[2,i] = float(rs.ring_count_mat[2,i])/float(length(rs.Dualnetwork.dual_node_list))
	end
	## mean
	rs.mu = 0.0
	for i in 1:size(rs.ring_stat_mat,2)
		x = float(rs.ring_stat_mat[1,i])
		y = float(rs.ring_stat_mat[2,i])
		rs.mu += y*x*1.0
	end
	rs.mu = rs.mu / 1.0 # area should be 1.0
	## standard deviation
	rs.sig = 0.0
	for i in 1:size(rs.ring_stat_mat,2)
		x = float(rs.ring_stat_mat[1,i])
		y = float(rs.ring_stat_mat[2,i])
		rs.sig += (x-rs.mu)*(x-rs.mu)*y*1.0
	end
	rs.sig = sqrt(rs.sig/1.0)
end

## calculate pi-matrix (for ring neighborhood statistics)
function calc_pi_matrix(rs::Ring_statistics)
	##
	rs.pi_matrix = zeros(length(rs.Dualnetwork.r_list),length(rs.Dualnetwork.r_list))
	cen_cntr = 1
	for cen_num in rs.r_list
		nei_cntr = 1
		for nei_num in rs.Dualnetwork.r_list
			## calculate entry of pi matrix
			sum_center_n_fold = 0
			pi_num = 0.0
			for center_node in rs.Dualnetwork.dual_node_list
				if center_node.type == cen_num
					sum_center_n_fold += 1
					sum_adj_nei = 0
					for adj_num in center_node.adjacent_node_list
						## add the neighbor to total sum
						if rs.Dualnetwork.dual_node_list[adj_num].type == nei_num
							sum_adj_nei += 1
						end
					end
					pi_num += float(sum_adj_nei) / float(length(center_node.adjacent_node_list))
				end
			end
			## save entry into pi matrix
			if sum_center_n_fold > 0
				rs.pi_matrix[cen_cntr,nei_cntr] = pi_num/float(sum_center_n_fold)
			end
			nei_cntr += 1
		end
		cen_cntr += 1
	end
end

## calculate Aboav-Weaire-Parameter
function calc_aboav_weaire(rs::Ring_statistics)
	rs.alpha_aboav_weaire_list = zeros(length(rs.r_list))
	for i in 1:length(rs.r_list)
		m = rs.r_list[i]
		mu_m = 0.0
		for ii in 1:length(rs.r_list)
			n = rs.r_list[ii]
			mu_m += rs.pi_matrix[i,ii]*n
		end
		if m != 6
			rs.alpha_aboav_weaire_list[i] = 1.0 - (m*mu_m-rs.mu*rs.mu-rs.sig*rs.sig) / (rs.mu*(m-rs.mu))
		end
	end
end

## objective function
function calc_chi(rs::Ring_statistics,num_c_ring=[4,5])
	a0_list = [2.0,2.0]
	a1 = 1.0
	## term for the ring neighborhood statistics
	rs.chi_alpha = 0.0
	for i in 1:length(num_c_ring)
		rs.chi_alpha += a0_list[i]*abs(rs.alpha_aboav_weaire_list[i] -
		                            rs.alpha_aboave_weaire_list_target[i])
	end
	## term for the ring statistics
	rs.chi_het = 0.0
	for n in rs.r_list
		pn = rs.ring_stat_mat[2,n]
		pt = rs.target_stat_vec[n]
		rs.chi_het += abs(pn - pt)/pt
	end
	rs.chi_het = a1*rs.chi_het
	## return total chi
	return rs.chi_alpha + rs.chi_het
end

## make decision using the Metropolis condition
function make_decision(rs::Ring_statistics,delta_chi,T=1.0e-4)
	# decision to accept a switch
	decision = true
	# probability of acceptance
	P = 0.0
	if delta_chi < 0.0
		P = 1.0
	else
		P = exp((-1.0)*delta_chi/T)
	end
	# create a random number between 0 and 1
	random_number = rand(Float64)
	if P > random_number
		decision = true
	else
		decision = false
	end
	return decision
end
