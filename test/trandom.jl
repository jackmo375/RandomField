include("../src/RandomField.jl")
using .RandomField

# choose parameters:
## remember, we need that 4|n_i for i=1,...,d
N = [400,400]		# number pixels per dimension
L = [20.0,80.0]		# real space box is [-l,l]^d
f(V::Vector{T} where T<:Real) = gaus_rf(V)
					# spectral density

# derived parameters
d = length(N) 		# dimension
H = (pi/2) .* N ./L
var_y(V::Vector{T} where T<:Real) = 2*(2pi)^(d/2)*prod(L)/pi^d*f(V)

# allocate memory
V = Cloud(d,N)
X = Cloud(d,N)
Y = Cscalar(V)
Z = Cscalar(X)

# initialize grids
init!(X,L)
init!(V,H)

draw_realization!(Y,var_y)

# Y ↦ Z
ift!(Z,Y)

# print:
tovtk(Y, "../data/v.y.vtk")
tocsv(slice(Y,1), "../data/v.y1.csv")
tovtk(Z, "../data/x.z.vtk")
tocsv(slice(Z,1), "../data/x.z1.csv")
