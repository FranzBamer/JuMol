## Read LRS results from file
#  empirical distribution density
######################################
######################################
######################################
# Fabric tensors
#
# (c) Franz Bamer, Ivy Wu September-2022
######################################


## calculate empirical distribution density
function calc_empirical_distribution_density(num_intersec, bond_list)
	dphi = 2*pi/float(num_intersec)
	hist = zeros(num_intersec)
	#
	for bond in bond_list
		# four quadrants
		alpha = atan(abs(bond[2]/bond[1]))
		if bond[1] >= 0 && bond[2] >= 0
			alpha += 0
		elseif bond[1] < 0 && bond[2] >= 0
			alpha = pi - alpha
		elseif bond[1] < 0 && bond[2] < 0
			alpha += pi
		else
			alpha = 2*pi - alpha
		end
		# find the angle for the histogram and add the counter
		for i in 1:num_intersec
			phi_low = dphi*(i-1)
			phi_high = dphi*i
			if alpha >= phi_low && alpha < phi_high
				hist[i] += 1
				break
			end
		end
	end
	# normalize histogram
	normalize_hist(hist, num_intersec)
	return hist
end


## extract fabric tensors of every local region
function calc_fabric_tensor_list(bond_list_list)
	# change list of matrices to list of list of vectors and calc fabric tensors
	F_mean_list = Vector()
	#lr.F_sig_list = Vector()
	for bond_list in bond_list_list
		F_mean = calc_fabric_tensor(bond_list)
		push!(F_mean_list, F_mean)
	end
	return F_mean_list
end

## extract fabric tensors of every local region componentwise standard deviation
function calc_fabric_tensor_list_sig(F_mean_list, bond_list_list)
	# change list of matrices to list of list of vectors and calc fabric tensors
	F_sig_list = Vector()
	cntr = 1
	for bond_list in bond_list_list
		F_sig = calc_fabric_tensor_sig(bond_list, F_mean_list[cntr])
		push!(F_sig_list, F_sig)
		cntr += 1
	end
	return F_sig_list
end

## calculate the mean of the fabric tensors in every local structure
function calc_fabric_tensor(bond_list)
	# mean of the fabric tensors
	F_mean = zeros(2,2)
	for bond_vec in bond_list
		F_mean += bond_vec*transpose(bond_vec)
	end
	F_mean = F_mean ./ float(length(bond_list))
	#
	return F_mean
end

## calculation of the standard deviation of the fabric tensor (component-wise)
function calc_fabric_tensor_sig(bond_list, F_mean)
	# standard deviations of the fabric tensors
	F_sig = zeros(2,2)
	for bond_vec in bond_list
		# calculate variance component-wise
		F_sig += (bond_vec*transpose(bond_vec) - F_mean) .^2
	end
	F_sig = F_sig ./ float(length(bond_list))
	# calc square root component-wise
	F_sig = F_sig.^(0.5)
	return F_sig
end
