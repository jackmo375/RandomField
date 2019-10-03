include("../src/RandomField.jl")
using .RandomField

# choose parameters:
## remember, we need that 4|n_i for i=1,...,d
N = [200,200]		# number pixels per dimension
L = [50.0,80.0]		# real space box is [-l,l]^d
f(V::Vector{T} where T<:Real) = rect_rf(V)

# derived parameters
d = length(N) 		# dimension
H = (pi/2) .* N ./L 

# allocate memory
V = Cloud(d,N)
X = Cloud(d,N)
F = Cscalar(V)
G = Cscalar(X)

# initialize grids
init!(X,L)
init!(V,H)

sample!(F,f)

ift!(G,F)

tocsv(slice(G,2), "../data/x.g1.csv")


