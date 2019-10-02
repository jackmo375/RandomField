function draw_realization!(
	result::Array{Complex{T}} where T<:AbstractFloat,
	g::Function,
	V::AbstractArray)

	@. result = randn(ComplexF64) * sqrt(g(V))
end

