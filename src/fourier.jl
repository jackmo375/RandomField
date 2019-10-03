"""
	test()

an function to check the module loads.
"""
function test()
	print("Hello world\n")
end


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
	if Y.cl.dim==1
		for i=1:Npoints(Y.cl)
			Y.values[i] *= (-1)^(i-1)
		end
	elseif Y.cl.dim==2
		for i=1:Npoints(Y.cl)
			Y.values[i] *= (-1)^(
				rem(i-1,Y.cl.N[1])
				+ div(i-1,Y.cl.N[1]))
		end
	elseif Y.cl.dim==3
		for i=1:Npoints(Y.cl)
			Y.values[i] *= (-1)^(
				rem(i-1,Y.cl.N[1])
				+ rem(div(i-1,Y.cl.N[1]),Y.cl.N[2])
				+ div(i-1,Y.cl.N[1]*Y.cl.N[2]))
		end
	else
		throw(ArgumentError("cloud dimension must be 1,2, or 3"))
	end
end

# ift: F ↦ G
function ift!(G::Cscalar, F::Cscalar)

	G.values .= F.values

	dims = Tuple(n for n in G.cl.N)

	phase_shift!(G)
	ifft!(reshape(G.values,dims))
	phase_shift!(G)

	for i in 1:G.cl.dim
		G.values .*= G.cl.N[i]/G.cl.points[end,i]
	end
	G.values .*= (pi/2)^(G.cl.dim/2)

end

# inverse Hankel?