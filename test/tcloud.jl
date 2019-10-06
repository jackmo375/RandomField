include("../src/RandomField.jl")
using .RandomField

# choose parameters:
## remember, we need that 4|n_i for i=1,...,d
N = [10,20,30]		# number pixels per dimension
L = [80.0,80.0,80.0]	# real space box is [-l,l]^d
f(V::Vector{T} where T<:Real) = rect_rf(V)

# derived parameters
d = dim(N) # dimension
H = (pi/2) .* N ./L 

# allocate memory
V = Cloud(d,N)
X = Cloud(d,N)
F = Cscalar(V)
G = Cscalar(X)

init!(X,L)
init!(V,H)

tovtk(X, "../data/x.cl.vtk")

Xx = slice(X,1)
tovtk(Xx, "../data/xx.cl.vtk")

Xy = slice(X,2)
tovtk(Xy, "../data/xy.cl.vtk")

Xz = slice(X,3)
tovtk(Xz, "../data/xz.cl.vtk")

sample!(F,f)

tovtk(F, "../data/v.F.vtk")

Fx = slice(F,1)
tovtk(Fx, "../data/v.f1.vtk")

Fy = slice(F,2)
tovtk(Fy, "../data/v.f2.vtk")

Fz = slice(F,3)
tovtk(Fz, "../data/v.f3.vtk")
tocsv(Fz, "../data/v.f3.csv")