#
#	Process
#
function create_grids(
	N::Int, 
	l::AbstractFloat)
	
	# compute derived parameters:
	Δt = 2*l/N
	h  = π*N/2/l
	Δv = 2*h/N

	# compute grids:
	t  = -l .+ collect(0:N-1).*Δt
	v  = -h .+ collect(0:N-1).*Δv

	return t, v
end

function phase_shift!(
	x::Vector{Complex{T}} where T<:AbstractFloat)

	for i in 1:length(x)
		x[i] *= (-1)^(i-1)
	end

	x
end

function ift!(
	G::Array{Complex{T}} where T<:AbstractFloat,
	F::Array{Complex{T}} where T<:AbstractFloat,
	N::Integer, 
	l::AbstractFloat)

	G .= F
	phase_shift!(ifft!(phase_shift!(G)))
	G .*= sqrt(pi/2)*(N/l)
end

#
#	Field
#

# Phase shift
function phase_shift!(Y::Cscalar)

	for i=1:Npoints(Y)
		Y.values[i] *= (-1)^sum(tovec(i,Y.cl.N))
	end
end

# ift: F ↦ G
function ift!(G::Cscalar, F::Cscalar)

	G.values .= F.values

	dims = Tuple(n for n in G.cl.N)

	phase_shift!(G)
	ifft!(reshape(G.values,dims))
	phase_shift!(G)

	L = G.cl.points[end,1:end]
	N = G.cl.N
	d = dim(G)
	G.values .*= (pi/2)^(d/2)*prod(N)/prod(L)

end

# inverse Hankel?