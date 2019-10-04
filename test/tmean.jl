include("../src/RandomField.jl")
using .RandomField

# choose parameters:
## remember, we need that 4|n_i for i=1,...,d
N = [400,400]		# number pixels per dimension
L = [20.0,80.0]		# real space box is [-l,l]^d
f(V::Vector{T} where T<:Real) = gaus_rf(V)
					# spectral density
n_real = 1000		# number of realizations 

# derived parameters
d = length(N) 		# dimension
H = (pi/2) .* N ./L
var_y(V::Vector{T} where T<:Real) = 2*(2pi)^(d/2)*prod(L)/pi^d*f(V)

# allocate memory
V = Cloud(d,N)
X = Cloud(d,N)
Y = Cscalar(V)
Z = Cscalar(X)
M_est = Cscalar(X)	# estimated mean

# initialize grids
init!(X,L)
init!(V,H)

# estimate mean
@. M_est.values = 0
for i=1:n_real
	draw_realization!(Y,var_y)
	ift!(Z,Y)	# Y ↦ Z
	@. M_est.values += Z.values
end
@. M_est.values /= n_real

# print:
tovtk(M_est, "../data/x.mu.vtk")
tocsv(slice(Z,1), "../data/x.z1.csv")
tocsv(slice(M_est,1), "../data/x.mu1.csv")
