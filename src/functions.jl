#
#	gaussian
#
# process:
gaus_rf(ν::T where T<:Real) = exp(-ν^2/2)
# field:
gaus_rf(V::Vector{T} where T<:Real) = gaus_rf(norm(V))


#
#	Rect / sinc
#
# process:
rect_rf(ν::T where T<:Real) = abs(ν) < 1 ? 1 : 0
sinc_rf(t::T where T<:Real) = (2/sqrt(2pi)) * (t!=0 ? sin(t)/t : 1)
# field:
rect_rf(V::Vector{T} where T<:Real) = rect_rf(norm(V))
sinc_rf(V::Vector{T} where T<:Real) = sinc_rf(norm(V))
