include("../src/RandomField.jl")
using .RandomField

# choose parameters:
## remember, we need that 4|n_i for i=1,...,d
N = [400,400]		# number pixels per dimension
L = [20.0,80.0]		# real space box is [-l,l]^d
f(V::Vector{T} where T<:Real) = gaus_rf(V)

# derived parameters
d = length(N) 		# dimension
H = (pi/2) .* N ./L 

# allocate memory
V = Cloud(d,N)
X = Cloud(d,N)
F = Cscalar(V)
G = Cscalar(X)
G_tru = Cscalar(X)

# initialize grids
init!(X,L)
init!(V,H)

sample!(F,f)

ift!(G,F)

sample!(G_tru,gaus_rf)

tocsv(slice(G,2), "../data/x.g2.csv")
tocsv(slice(G_tru,2), "../data/x.gtrue2.csv")


