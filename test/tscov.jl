include("../src/RandomField.jl")
using .RandomField

# choose parameters:
## remember, we need that 4|n_i for i=1,...,d
N = [400,400]		# number pixels per dimension
L = [20.0,80.0]		# real space box is [-l,l]^d
f(V::Vector{T} where T<:Real) = gaus_rf(V)
					# spectral density
n_real = 100		# number of realizations 

# derived parameters
d = length(N) 		# dimension
H = (pi/2) .* N ./L
var_y(V::Vector{T} where T<:Real) = 2*(2pi)^(d/2)*prod(L)/pi^d*f(V)

# allocate memory
V = Cloud(d,N)
X = Cloud(d,N)
Y = Cscalar(V)
Z = Cscalar(X)
C_est = Cscalar(X)	# estimated stationary covariance
C_tru = Cscalar(X)	# true stationary covariance

# initialize grids
init!(X,L)
init!(V,H)

# estimate stationary covariance
## (we assume here that the mean is 0 to save time)
@. C_est.values = 0
i_o = get_originindex(X)
for i=1:n_real
	draw_realization!(Y,var_y)
	ift!(Z,Y)	# Y ↦ Z
	@. C_est.values += 0.5 * Z.values * conj(Z.values[i_o])
end
@. C_est.values /= n_real

sample!(C_tru,gaus_rf)

# print:
tovtk(C_est, "../data/x.c.vtk")
#tocsv(slice(Z,1), "../data/x.z1.csv")
tocsv(slice(C_est,1), "../data/x.c1.csv")
tocsv(slice(C_tru,1), "../data/x.ctru1.csv")