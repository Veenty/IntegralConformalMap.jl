
using LinearOperators
using Krylov
using FFTW

export ConformalMapConstructor


function ConformalMapConstructor( γ; L=1. ,solvertype = :gmres ,  rtol = 1e-13, atol = 1e-13 , α = 5)

    γ′ = DerivativeFourier(γ, L )
    γ″ = Derivative2Fourier(γ, L )

    if α ≈ 1

        return ConformalMapConstructor(γ, γ′, γ″; solvertype = solvertype ,  rtol = rtol, atol = atol)


    end

    γ_over = OversamplingFourier(γ, α  )
    γ′_over = OversamplingFourier(γ′, α )
    γ″_over = OversamplingFourier(γ″, α)



    return ConformalMapConstructor(γ_over, γ′_over, γ″_over; solvertype = solvertype ,  rtol = rtol, atol = atol)

end


function ConformalMapConstructor(γ, γ′, γ″; solvertype = :gmres ,  rtol = 1e-13, atol = 1e-13)

    #Assume that N already has the length of the discretization wanted to evaluate the conformal map 

    #The associated density of the integral equation could be of lower bandwidth than 
    #the one of conformal map. I hope the nice parametrization that I chosed 

    N = length(γ)

    v = abs.(γ′)
    κ = (γ′ .⊗ γ″) ./ v.^3
    n̂ = -1im*γ′ ./ v

    w = ones(N)/N
    ẇ = w .* abs.( γ′)

    b = -log.(abs.(γ))

    
    if solvertype == :gmres

        D_op = LinearOperator(Float64, N, N, false, false,
        (res, φ, α, β) -> D_operator!( res ,φ ,γ, ẇ , κ, n̂, α, β ), 
                    nothing,       
                    nothing)

        Problem = -(1/2)*opEye(N,N) + D_op

        (φ, stats) = gmres( Problem, b, b; rtol = rtol, atol = atol)


    elseif solvertype == :direct 

        Dφ = similar(φ)
        Dop = φ -> D_operator( φ ,γ_disc, ẇ , κ_disc,  n̂_disc )
        MD = MatrixFormation(N,N, Dop, Float64)

        φ = (-I*1/2 + MD)\b

    end

    #Convergence Criteria 
    # φ̂ = fft(φ)

    h = zeros(ComplexF64, N)
    Generate_holomorphic_data(h, φ, 1., γ, w, γ′  )

    #Here I might want to do oversampling of h, γ and γ′
    γ̃ =  γ.* exp.(h)
    γ̃ₛ = γ′.* exp.(h) + DerivativeFourier(h, 1  ).* γ̃

    f = (F,ζ) -> C_operator!(F, ζ, γ ,γ̃, γ̃ₛ , w)
    f′ = (F,ζ) -> C_operator!(F, ζ, γ′./ γ̃ₛ ,γ̃, γ̃ₛ , w)
    f⁻¹ = (F,z) -> (C_operator!(F, z, h ,γ, γ̃ₛ , w), z.*exp.(F))[2]


    return f, f′, f⁻¹

end