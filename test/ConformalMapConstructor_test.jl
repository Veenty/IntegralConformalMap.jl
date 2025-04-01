
using IntegralConformalMap
using Revise
using LinearAlgebra
using CairoMakie
using ComplexPlots
using LinearOperators
# using Krylov 
using FFTW
using BenchmarkTools

ϵ = 0.5
k = 5

γ = s-> exp.(2π*1im*s)*(1 +ϵ*cos(2π*k*s)  )
# γₛ = s-> 2π*1im*exp.(2π*1im*s)*(1 + ϵ*cos(2π*k*s)  + 1im*ϵ*k*sin(2π*k*s)  )
# γₛₛ = s-> -4*π^2*exp.(2π*1im*s)*(1 + ϵ*(k^2 +1)*cos(2π*k*s)  + 2im*ϵ*k*sin(2π*k*s)  )


N = 100

γ_disc = γ.(range(0,1,length=N+1)[1:end-1])

plot( fft(γ_disc) .|> log10∘abs) 

f, f′, f⁻¹ =  ConformalMapConstructor( γ_disc; L=1. , α = 1)

# @benchmark ConformalMapConstructor( $γ_disc; L=1. , α = 1)


ρ = sqrt.(collect(0:0.1:1)[2:end])
ϕ = collect(0:0.001:2π)[1:end-1]
ζ = ρ*exp.(2π*1im*ϕ')
F = similar(ζ )


# @benchmark f($F,$ζ)

scatter(vec(F), markersize = 4, color = :blue)
