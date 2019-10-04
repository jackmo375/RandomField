include("../src/RandomField.jl")
using .RandomField

# chose input parameters:
n 		= 1000		# number of grid points
n_real	= 10000		# number of realizations
l 		= 30.0		# real space box is [-l,l]
f 		= gaus_rf	# spectral density

var_y(ν) = 2*sqrt(2pi)*(l/pi)*f(ν)

# create structures (allocate memory)
Y = Vector{ComplexF64}(undef,n)	# fourier realization
Z = Vector{ComplexF64}(undef,n)	# real realization

M_est = zeros(ComplexF64,n) # est. mean
C_est = zeros(ComplexF64,n)	# est. stationary covariance
F_est = zeros(ComplexF64,n)	# est. spectral density
C_tru = Vector{ComplexF64}(undef,n) # true stationary covariance
F_tru = Vector{ComplexF64}(undef,n)	# true spectral density

T, V = create_grids(n,l)

# draw realization
draw_realization!(Y,var_y,V)
ift!(Z,Y,n,l)

# zeroth order moment
for i=1:n_real
	draw_realization!(Y,var_y,V)
	ift!(Z,Y,n,l)
	M_est .+= Z
end
M_est ./= n_real

# first order moments
for i=1:n_real
	draw_realization!(Y,var_y,V)
	ift!(Z,Y,n,l)
	@. F_est += sqrt(1/8pi)*(pi/l) * Y * conj(Y)
	@. C_est += 0.5 * Z * conj(Z[div(n,2)])
end
F_est ./= n_real
C_est ./= n_real

@. F_tru = f(V)
ift!(C_tru,F_tru,n,l)

# print to file:
out = open("../data/output.csv", "w")
print(out, "T, V, ReY, ReZ, ImZ, M_est, C_est, C_tru, F_est, F_tru\n")
for i in 1:n
	print(out,
		"$(T[i]),"*
		"$(V[i]),"*
		"$(real(Y[i])),"*
		"$(real(Z[i])),"*
		"$(imag(Z[i])),"*
		"$(real(M_est[i])),"*
		"$(real(C_est[i])),"*
		"$(real(C_tru[i])),"*
		"$(real(F_est[i])),"*
		"$(real(F_tru[i]))\n")
end
close(out)