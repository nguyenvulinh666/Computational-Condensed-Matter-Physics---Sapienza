import Random as rd
import Printf
N::Int64 = 2                            # Number of particles (2 = Helium)
R = Array{Float64,2}(undef,3,N)         # Position of the electron(s)
randR = Array{Float64,2}(undef,3,N)     # Random change of position
delta::Float64 = 1                      # Size of the step 
RP = Array{Float64,2}(undef,3,N)        # Proposed changed position of the electron(s)
p::Float64 = 0                          # Probability of making the proposed change
k::Int64 = 1                            # Number of variational parameters
β::Vector{Float64} = [1.3]              # Variational parameters
# Function that returns the probability distribution at a certain point for certain parameters
function Ψ²(P::Matrix{Float64},α::Vector{Float64})::Float64
    return Base.Math.exp(-2*α[1].*sum(sqrt.(sum(P.*P,dims=1))))
end
function fEP(P::Matrix{Float64},α::Vector{Float64})::Float64
    return sum(-2.0./sqrt.(sum(P.*P,dims=1))) + 1/sqrt.(sum((P[:,2]-P[:,1]).^2))
end
function fEK(P::Matrix{Float64},α::Vector{Float64})::Float64
    return -0.5*(sum(α).^2-2*α[1]./sqrt(sum(P[:,1].^2))) -0.5*(α[1].^2-2*α[1]./sqrt(sum(P[:,2].^2)))
end
EK::Float64 = 0                         # Kinetic Energy
EP::Float64 = 0                         # Potential Energy
ET::Float64 = 0                         # Total Energy
print("i,beta,EK,EP,ET,p")
for l in 1:N
    Printf.@printf(",r%d",l)
    Printf.@printf(",x%d,y%d,z%d",l,l,l)
end
Printf.@printf("\n")
for beta in 1.6:0.05:1.8
    global β = [beta]
    rd.rand!(R); global R *= 10             # Initial positions
    n::Int64 = 1000000                       # Number of mc steps
    for i in 1:n
        rd.rand!(randR)                     # Generate random change in the positions
        global RP = R+delta.*(randR.-0.5)   # Obtain the proposed change
        global p = min(1,Ψ²(RP,β)/Ψ²(R,β))  # Probability for the proposed change
        if p > rand()
            global R  = copy(RP)
            global EK = fEK(R,β)
            global EP = fEP(R,β)
            global ET = EK+EP
        end
        Printf.@printf("%d,%f,%f,%f,%f,%f",i,sum(beta),EK,EP,ET,p)
        for l in 1:N
            Printf.@printf(",%f",sqrt.(sum(R.*R,dims=2))[l])
            Printf.@printf(",%f,%f,%f",R[1,l],R[2],R[3,l])
        end
        Printf.@printf("\n")
    end   
end
