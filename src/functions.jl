#
#	gaussian
#
gaus_rf(ν) = exp(-ν^2/2)


#
#	Rect / sinc
#
rect_rf(ν) = abs(ν) < 1 ? 1 : 0
sinc_rf(t) = (2/sqrt(2pi)) * (t!=0 ? sin(t)/t : 1)
