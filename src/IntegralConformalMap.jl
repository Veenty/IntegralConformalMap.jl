module IntegralConformalMap

include("IntegralEquation.jl")
include("ConformalMapConstructor.jl")

# Write your package code here.

export MatrixFormation, D_operator!, D_operator, Generate_holomorphic_data, C_operator!, DerivativeFourier
export C_operator_direct
export ⊙, ⊗
export ConformalMapConstructor

end
