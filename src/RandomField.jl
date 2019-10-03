module RandomField

using FFTW, Random, LinearAlgebra

export
	# from functions.jl
	gaus_rf,
	rect_rf,
	sinc_rf,

	# from cloud.jl:
	Cloud,
	init!,
	tovtk,
	tocsv,
	Cscalar,
	sample!,
	slice,

	# from fourier.jl:
	test,
	create_grids,
	#phase_shift,
	ift!,

	# from process.jl
	draw_realization!

	# from FFTW
	#ifft

include("functions.jl")

include("cloud.jl")

include("fourier.jl")

include("process.jl")

end # of RandomField