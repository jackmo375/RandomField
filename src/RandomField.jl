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
	dim,
	surface,

	# from fourier.jl:
	create_grids,
	ift!,

	# from random.jl
	draw_realization!

include("functions.jl")

include("cloud.jl")

include("fourier.jl")

include("random.jl")

end # of RandomField
