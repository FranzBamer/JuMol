## read files
######################################
######################################
######################################


using DelimitedFiles



mutable struct Read
	lines::Vector
	line_start_list::Vector
	Read() = new()
	## write file
	output_file
	## translate when importing from lammps
	x_translate::Float64
	y_translate::Float64
	z_translate::Float64
end





##############################################
#### reading the time history of a lammpsrj
##############################################
function get_lines(reader::Read, filename)
	## define lines and line_start_list
	reader.lines = Vector{String}()
	reader.line_start_list = Vector{Int64}()
	## read lines
	file = open(filename,"r")
	reader.lines = readlines(file)
	close(file)
	## find line start list of every timestep
	cntr = 1
	while cntr < length(reader.lines)
		if reader.lines[cntr] == "ITEM: TIMESTEP"
			push!(reader.line_start_list,cntr)
		end
		cntr += 1
	end
end

#### get the time step
function read_lammpstrj(reader::Read, structure, line_start, triclinic=false)
	#### begin to read the file
	## description of the time step
	cntr = line_start + 1
	## timestep
	num_timestep = parse(Int,reader.lines[cntr],base=10)
	## description of number of atoms
	cntr += 2
	## number of atoms
	structure.noa = parse(Int,reader.lines[cntr],base=10)
	## box boundaries
	cntr += 2
	line_vec = readdlm(IOBuffer(reader.lines[cntr]))
	reader.x_translate = line_vec[1]
	structure.box.lx = line_vec[2] - reader.x_translate
	if triclinic
		structure.box.lxy = line_vec[3]
	else
		structure.box.lxy = 0.0
	end
	cntr += 1
	line_vec = readdlm(IOBuffer(reader.lines[cntr]))
	reader.y_translate = line_vec[1]
	structure.box.ly = line_vec[2] - reader.y_translate
	if triclinic
		structure.box.lxz = line_vec[3]
	else
		structure.box.lxz = 0.0
	end
	cntr += 1
	line_vec = readdlm(IOBuffer(reader.lines[cntr]))
	reader.z_translate = line_vec[1]
	structure.box.lz = line_vec[2] - reader.z_translate
	if triclinic
		structure.box.lyz = line_vec[3]
	else
		structure.box.lyz = 0.0
	end
	#### reading all atoms
	cntr += 1
	for i in 1:structure.noa
		line_vec = readdlm(IOBuffer(reader.lines[cntr+i]))
		num_atom = Int(line_vec[1])
		type_atom = Int(line_vec[2])
		push!(structure.atom_list,Atom(num_atom,type_atom))
		structure.atom_list[i].pos = [Float64(line_vec[3])-reader.x_translate,Float64(line_vec[4])-reader.y_translate,Float64(line_vec[5])-reader.z_translate]
		structure.atom_list[i].vel = [Float64(line_vec[6]),Float64(line_vec[7]),Float64(line_vec[8])]
		structure.atom_list[i].vel = [0.0,0.0,0.0]
	end
end

# ############################################
# #### get the time step from Spencer
# function read_lammpstrj_spencer(reader::Read, structure, line_start)
# 	#### begin to read the file
# 	## description of the time step
# 	cntr = line_start + 1
# 	## timestep
# 	num_timestep = parse(Int,reader.lines[cntr],base=10)
# 	## description of number of atoms
# 	cntr += 2
# 	## number of atoms
# 	structure.noa = parse(Int,reader.lines[cntr],base=10)
# 	## box boundaries
# 	cntr += 2
# 	line_vec = readdlm(IOBuffer(reader.lines[cntr]))
# 	reader.x_translate = line_vec[1]
# 	structure.box.lx = line_vec[2] - reader.x_translate
# 	cntr += 1
# 	line_vec = readdlm(IOBuffer(reader.lines[cntr]))
# 	reader.y_translate = line_vec[1]
# 	structure.box.ly = line_vec[2] - reader.y_translate
# 	cntr += 1
# 	structure.box.lz = 10.0
# 	#### reading all atoms
# 	cntr += 1
# 	for i in 1:structure.noa
# 		line_vec = readdlm(IOBuffer(reader.lines[cntr+i]))
# 		num_atom = Int(line_vec[1])
# 		type_atom = Int(line_vec[2])
# 		push!(structure.atom_list, Atom(i, type_atom))
# 		structure.atom_list[i].pos = [Float64(line_vec[3])-reader.x_translate, Float64(line_vec[4])-reader.y_translate, 0.0]
# 		#structure.atom_list[i].vel = [Float64(line_vec[6]),Float64(line_vec[7]),Float64(line_vec[8])]
# 		structure.atom_list[i].vel = [0.0,0.0,0.0]
# 	end
# end

##############################################
#### load the box from a box input file
##############################################
function read_box_lammps(structure,filename)
	## read all the lines from the file
	file = open(filename,"r")
	lines = readlines(file)
	close(file)

	## Read number of atoms
	cntr =3
	line_vec = readdlm(IOBuffer(lines[cntr]))
	structure.noa = line_vec[1]

	#### reading all atoms
	cntr = 17
	for i in 1:structure.noa
		line_vec = readdlm(IOBuffer(lines[cntr+i]))
		num_atom = Int(line_vec[1])
		type_atom = Int(line_vec[2])
		push!(structure.atom_list,Atom(num_atom,type_atom))
		structure.atom_list[i].pos = line_vec[3:5]
	end
	println(structure.atom_list)
end

## open a file for output
# write_choice: w ... write, a ... append, r ... read
function open_file(reader::Read, filename::String, write_choice::String)
	reader.output_file = open(filename, write_choice)
end

## write a lammpstj file
function write_box(reader::Read,structure,num_step)
	write(reader.output_file,"ITEM: TIMESTEP\n")
	write(reader.output_file,string(num_step,"\n"))
	write(reader.output_file,"ITEM: NUMBER OF ATOMS\n")
	write(reader.output_file,string(structure.noa,"\n"))
	write(reader.output_file,"ITEM: BOX BOUNDS xy xz yz pp pp pp\n")
	write(reader.output_file,string(0.0," ", structure.box.lx, " ", structure.box.lxy,"\n"))
	write(reader.output_file,string(0.0," ", structure.box.ly, " ", structure.box.lyz,"\n"))
	write(reader.output_file,string(0.0," ", structure.box.lz, " ", structure.box.lxz,"\n"))
	write(reader.output_file,"ITEM: ATOMS id type x y z vx vy vz\n")
	## write atom positions and velocities
	for i in 1:structure.noa
		atom = structure.atom_list[i]
		line = string(atom.number," ",atom.type," ",atom.pos[1]," ",atom.pos[2], " ",atom.pos[3],
		              " ", atom.vel[1], " ", atom.vel[2], " ", atom.vel[3],"\n")
		write(reader.output_file,line)
	end
end

## write a twodsilica file
function write_lmp_box(reader::Read, structure, num_step)
	write(reader.output_file,"#input file\n")
	write(reader.output_file,"\n")
	write(reader.output_file,string(structure.noa)," atoms\n")
	write(reader.output_file,"2 atom types\n")
	write(reader.output_file,"\n")
	write(reader.output_file,string(0.0)," ", string(structure.box.lx), " xlo xhi\n")
	write(reader.output_file,string(0.0)," ", string(structure.box.ly), " ylo yhi\n")
	write(reader.output_file,"-10.0 10.0 zlo zhi\n")
	write(reader.output_file,string(structure.box.lxy)," 0.0 0.0 xy xz yz\n") #FIXME
	write(reader.output_file,"\n")
	write(reader.output_file,"Masses\n")
	write(reader.output_file,"\n")
	write(reader.output_file,"1 28.0855 #Si\n")
	write(reader.output_file,"2 15.9994 #O\n")
	write(reader.output_file,"\n")
	write(reader.output_file,"Atoms\n")
	write(reader.output_file,"\n")
	## write atom positions and velocities
	for i in 1:structure.noa
		atom = structure.atom_list[i]
		line = string(atom.number," ",atom.type," ",atom.pos[1]," ",atom.pos[2], " 0.0\n")
		write(reader.output_file,line)
	end
end

## write a twodsilica group file
function write_lmp_box_group(reader::Read, structure, num_step::Int64, num_group::Int64)
	# number of group atoms
	cntr_group = 0
	for atom in structure.atom_list
		if atom.group == num_group
			cntr_group += 1
		end
	end
	# write header
	write(reader.output_file,"#input file\n")
	write(reader.output_file,"\n")
	write(reader.output_file,string(cntr_group)," atoms\n")
	write(reader.output_file,"2 atom types\n")
	write(reader.output_file,"\n")
	write(reader.output_file,string(0.0)," ", string(structure.box.lx), " xlo xhi\n")
	write(reader.output_file,string(0.0)," ", string(structure.box.ly), " ylo yhi\n")
	write(reader.output_file,"-10.0 10.0 zlo zhi\n")
	write(reader.output_file,string(structure.box.lxy)," 0.0 0.0 xy xz yz\n") #FIXME
	write(reader.output_file,"\n")
	write(reader.output_file,"Masses\n")
	write(reader.output_file,"\n")
	write(reader.output_file,"1 28.0855 #Si\n")
	write(reader.output_file,"2 15.9994 #O\n")
	write(reader.output_file,"\n")
	write(reader.output_file,"Atoms\n")
	write(reader.output_file,"\n")
	## write atom positions and velocities
	for atom in structure.atom_list
		if atom.group == num_group
			line = string(atom.number," ",atom.type," ",atom.pos[1]," ",atom.pos[2], " 0.0\n")
			write(reader.output_file,line)
		end
	end

end

## close file for output
function close_file(reader::Read)
	close(reader.output_file)
end
