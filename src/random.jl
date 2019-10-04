#
#	Process
#

function draw_realization!(
	result::Array{Complex{T}} where T<:AbstractFloat,
	g::Function,
	V::AbstractArray)

	@. result = randn(ComplexF64) * sqrt(g(V))
end

#
#	Field
#
function draw_realization!(Y::Cscalar,var::Function)
	for i=1:Npoints(Y)
		Y.values[i] = randn(ComplexF64) * sqrt(var(Y.cl.points[i,1:end]))
	end
end

