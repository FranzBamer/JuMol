## Read LRS results from file
#  general functions
######################################
######################################
######################################
#
#
# (c) Franz Bamer, Ivy Wu September-2022
######################################





## function matrix to vector
function mat2vec_list(mat)
	vec_list = Vector()
	for i in 1:size(mat,2)
		vec = mat[:,i]
		push!(vec_list, vec)
	end
	return vec_list
end

## function to normalize a list of vectors
function normalize_list(bond_list)
	n_bond_list = Vector()
	l_bond_list = Vector()
	for i in 1:length(bond_list)
		l_bond = norm(bond_list[i])
		n_bond = bond_list[i] ./ norm(bond_list[i])
		push!(l_bond_list, l_bond)
		push!(n_bond_list, n_bond)
	end
	return n_bond_list, l_bond_list
end

## function to normalize the histogram
function normalize_hist(hist, num_seg)
	sum = 0.0
	for val in hist
		sum += val
	end
	#
	for i in 1:length(hist)
		hist[i] = hist[i] / (sum*2pi/num_seg)
	end
end

## calculate the eigensolution of a 2x2 matrix
function calc_eigensolution_two_times_two(F::Array)
	lambdavec = zeros(2)
	sqroot = sqrt((F[1,1]-F[2,2])^2 * 0.25 + F[1,2]*F[2,1])
	lambdavec[1] = (F[1,1]+F[2,2])*0.5 + sqroot
	lambdavec[2] = (F[1,1]+F[2,2])*0.5 - sqroot
	vec1 = zeros(2)
	vec1[1] = 1.0
	vec1[2] = -F[2,1] / (F[2,2]-lambdavec[1])
	vec2 = zeros(2)
	vec2[1] = 1.0
	vec2[2] = -F[2,1] / (F[2,2]-lambdavec[2])
	vec1 = vec1 ./ norm(vec1)
	vec2 = vec2 ./ norm(vec2)
	eigen_direc_mat = zeros(2,2)
	eigen_direc_mat[:,1] = vec1
	eigen_direc_mat[:,2] = vec2
	#
	return lambdavec, eigen_direc_mat
end
