include("../src/RandomField.jl")
using .RandomField

# chose input parameters:
n    = 1000			# number of grid points
l    = 80.0			# real space box is [-l,l]
f = rect_rf	# spectral density

# create structures
F = Vector{ComplexF64}(undef,n)
G = Vector{ComplexF64}(undef,n)

T, V = create_grids(n, l)

F .= f.(V)

ift!(G,F,n,l)

# print to file:
out = open("../data/output.csv", "w")
print(out, "T, V, ReF, ReG, expected\n")
for i in 1:n
	print(out,
		"$(T[i]),"*
		"$(V[i]),"*
		"$(real(F[i])),"*
		"$(real(G[i])),"*
		"$(sinc_rf(T[i]))\n")
end
close(out)