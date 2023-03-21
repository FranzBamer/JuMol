## module defintion of JuMol
######################################
######################################
######################################


module Jumol
    include("../src/sys/sys.jl")
    include("../src/structure/Structure.jl")
    include("../src/vis/Vis2d.jl")
    include("../src/calc/min.jl")
    include("../src/calc/int.jl")
    include("../src/calc/calc.jl")
    include("../src/calc/units.jl")
    include("../src/structure/Modifier.jl")
    include("../src/deform/Deform.jl")
    include("../src/strategies/Box_relax.jl")
    include("../src/calc/hessian.jl")
    include("../src/strategies/Gen_lattices.jl")
    include("../src/strategies/dual_structure/Dual_network.jl")
    include("../src/strategies/dual_structure/Dual_modifier.jl")
    include("../src/strategies/phys_structure/Phys_network.jl")
    include("../src/vis/Network_vis.jl")
    include("../src/strategies/Network.jl")
    include("../src/deform/Disp_field.jl")
    include("../src/deform/Strain_tensor_field.jl")
    include("../src/strategies/LYS/LYS.jl")
    include("../src/strategies/LYS/LYS_Read.jl")
    include("../src/strategies/LRS/LRS.jl")
    include("../src/strategies/LRS/LRS_Read.jl")
end
