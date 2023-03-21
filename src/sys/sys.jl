## system functions
######################################
######################################
######################################
#
# (c) Franz Bamer, Mai-2022
######################################


## set the current directory for calculations
function set_the_current_path()
	path_to_current_file = Base.source_path()
	del_char = "."
	path = collect(path_to_current_file)
	while path[end] != '/'
		char_del = path[end]
		num_delete = length(path)
		path = deleteat!(path,num_delete)
	end
	#path_to_current_file = "/Users/franzbamer/Library/Mobile Documents/com~apple~CloudDocs/evaluations_cloud/jumol_eval/2022_proj_LYS_40x48_scan/"
	path_to_current_file = String(path)
	cd(path_to_current_file)
	println("switching to")
	run(`pwd`)
end
