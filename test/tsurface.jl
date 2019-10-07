include("../src/RandomField.jl")
using .RandomField

# choose parameters:
## remember, we need that 4|n_i for i=1,...,d
N = [100,100]		# number pixels per dimension
L = [5.0,5.0]	# real space box is [-l,l]^d
f(X::Vector{T} where T<:Real) = 10*gaus_rf(X)

# derived parameters
d = dim(N) # dimension
H = (pi/2) .* N ./L 

# allocate memory
X = Cloud(d,N)
G = Cscalar(X)

init!(X,L)

sample!(G,f)

tovtk(surface(G), "../data/surface.vtk")