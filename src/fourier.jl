"""
	test()

an function to check the module loads.
"""
function test()
	print("Hello world\n")
end

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
	x::AbstractArray)

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