
# using LinearAlgebra
using FFTW

function ⊙(x,y)

    return real(x)*real(y) + imag(x)*imag(y)

end

function ⊗(x,y)

    return real(x)*imag(y) - imag(x)*real(y)

end

function MatrixFormation(N,M, func, T)

    A = zeros(T, N, M)
    eⱼ = zeros(T, M)

    for j ∈ 1:M

        eⱼ[j]= 1

        A[:,j] .= func(eⱼ)

        eⱼ[j] = 0

    end

    return A



end




function D_operator!( Dφ ,φ ,γ, ẇ , κ, n̂ )

    Dφ .= 0

    for j ∈ eachindex(φ)
        
        for i ∈ eachindex(γ)

            Δx = γ[i] - γ[j]

            if abs(Δx) ≈ 0

                Dφ[i] -= (1/2π) * κ[i]/2*φ[j]*ẇ[j]

            else

                Dφ[i] -= -(1/2π) * Δx ⊙ n̂[j]  / abs(Δx)^2 *ẇ[j]*φ[j]

            end
                 
            
        end

    end




end

function D_operator!( Dφ ,φ ,γ, ẇ , κ, n̂, α, β )

    if β == 0
        Dφ .= 0

    else

        Dφ .=β*Dφ
   
    end

     

    for j ∈ eachindex(φ)
        
        for i ∈ eachindex(γ)

            Δx = γ[i] - γ[j]

            if abs(Δx) ≈ 0

                Dφ[i] -= (α/2π) * κ[i]/2*φ[i]*ẇ[i]

            else

                Dφ[i] -= -(α/2π) * Δx ⊙ n̂[j]  / abs(Δx)^2 *ẇ[j]*φ[j]

            end
                 
            
        end

    end




end


function D_operator( φ ,γ, ẇ , κ, n̂ )

    Dφ = similar(φ)

    D_operator!( Dφ ,φ ,γ, ẇ , κ, n̂ )

    return Dφ
end


function C_operator_direct(x, φ ,γ, γ′ , w)

    Cφ = 0.0im

    for i ∈ eachindex(γ)
        Δz = γ[i] - x


        Cφ += φ[i]*w[i]*γ′[i]/Δz

    end

    return Cφ/(2π*1im)

end



function C_operator!(Cf, X , f ,γ, γ′ , w)

    #Cauchy Operator for holomorphic functions
    #We use barycentric Approximation  

    Cf .= 0im 

    for j ∈ eachindex(X)
        Int1 = 0
        Int2 = 0
        for i ∈ eachindex(γ )
            Δz = γ[i] - X[j] 
            if abs(Δz) ≈ 0
                Int1 = f[j]
                Int2 = 1
                break
            end

            Int1 += f[i]*w[i]*γ′[i]/Δz
            Int2 += w[i]*γ′[i]/Δz
    
        end
        Cf[j] = Int1/Int2

    end


end



function DerivativeFourier(τ, L )

    N = length(τ)
    τ̂ = fft(τ)
    k = zeros(N)
    k[1:N÷2+1] .=  collect(0:N÷2)
    k[N÷2+2:end] .=   collect((N÷2+1:(N-1))).-(N)
    k[N÷2+1] = k[N÷2+1]*(N%2) 

    # display(k)

    return 2π/L*ifft(1im*(k.*τ̂) )


end


function Derivative2Fourier(τ, L )

    N = length(τ)
    τ̂ = fft(τ)
    k = zeros(N)
    k[1:N÷2+1] .=  collect(0:N÷2)
    k[N÷2+2:end] .=   collect((N÷2+1:(N-1))).-(N)
    # k[N÷2+1] = k[N÷2+1]*(N%2) 

    # display(k)

    return -(2π/L)^2*ifft((k.^2 .*τ̂) )


end 


function OversamplingFourier(τ, α  )

    N = length(τ)
    τ̂ = fft(τ)
    M = floor(N*α)
    τ̂_padded = zeros(eltype(τ̂), M)
    τ̂_padded[1:N÷2+1] = τ̂[1:N÷2+1]
    τ̂_padded[N÷2+1] = τ̂_padded[N÷2+1]*(1/2 + (N%2)/2)
    tail_size = (N÷2 -1)*(1 - (N%2)) + N÷2*(N%2)
    τ̂_padded[end-tail_size+1:end] = τ̂[N÷2+2:end]
    τ̂_padded[end-tail_size] = τ̂_padded[N÷2+1]*(1 - (N%2)) 

    return ifft(τ̂_padded)*(M/N)
    


end



function Generate_holomorphic_data(f, τ, L, γ ,w, γ′  )

    τ′ = DerivativeFourier(τ, L)

    for k ∈ eachindex(f)
        f[k] = -τ[k]

        for j ∈ eachindex(τ)
            Δτ = τ[j] - τ[k]
            Δz = γ[j] - γ[k]

            if k == j
                
                f[k] -= 1/(2π*1im) * τ′[k]*w[k]
            else

                f[k] -= 1/(2π*1im) * Δτ/Δz *w[j]*γ′[j]

            end

        end
        # f[i] = -f[i]


    end

end


