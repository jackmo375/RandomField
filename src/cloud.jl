#
#	Types
#
mutable struct Cloud
	dim::Int64
	N::Vector{Int64}
	points::Array{Float64,2}
end

mutable struct SubCloud
	cl::Cloud
	N::Vector{Int64}
	verts::Vector{Int64}
end

# complex scalar field
mutable struct Cscalar
	cl::Cloud
	values::Vector{ComplexF64}
end

mutable struct SubCscalar
	cs::Cscalar
	N::Vector{Int64}
	verts::Vector{Int64}
end

#
#	Functions
#

function getverts(cl::Cloud, axis::Integer)
	if cl.dim==1
		i_start = 1
		i_stop  = cl.N[1]
		verts   = collect(i_start:i_stop)
		N 		= cl.N
	elseif cl.dim==2
		if axis==1
			i_start = 1 + div(cl.N[2]*cl.N[1],2)
			i_stop  = i_start + cl.N[1] - 1
			verts   = collect(i_start:i_stop)
			N = [cl.N[1],1]
		elseif axis==2
			i_start = 1 + div(cl.N[1],2)
			i_stop  = cl.N[1]*cl.N[2] - 1
			verts   = collect(i_start:cl.N[1]:i_stop)
			N = [1,cl.N[2]]
		else
			throw(ArgumentError("your cloud has < 3 dimensions"))
		end
	elseif cl.dim==3
		if axis==1
			i_start = 1 + div(cl.N[3]*cl.N[2]*cl.N[1] + cl.N[2]*cl.N[1],2)
			i_stop  = i_start + cl.N[1] - 1
			verts   = collect(i_start:i_stop)
			N = [cl.N[1],1,1]
		elseif axis==2
			i_start = 1 + div(cl.N[3]*cl.N[2]*cl.N[1] + cl.N[1],2)
			i_stop  = i_start + cl.N[2]*cl.N[1] - 1
			verts   = collect(i_start:cl.N[1]:i_stop)
			N = [1,cl.N[2],1]
		elseif axis==3
			i_start = 1 + div(cl.N[2]*cl.N[1] + cl.N[1],2)
			i_stop  = i_start + cl.N[3]*cl.N[2]*cl.N[1] - 1
			verts   = collect(i_start:cl.N[2]*cl.N[1]:i_stop)
			N = [1,1,cl.N[3]]
		else
			throw(ArgumentError("you can only slice in one of 3 dimensions"))
		end
	else
		throw(ArgumentError("you can only slice 2 and 3 dimension clouds"))
	end

	return N, verts
end

function slice(cl::Cloud, axis::Integer)
	
	N, verts = getverts(cl,axis)

	SubCloud(cl,N,verts)
end

function slice(cs::Cscalar, axis::Integer)

	N, verts = getverts(cs.cl, axis)

	SubCscalar(cs, N, verts)
end

# Constructors
function Cloud(
	d::Integer,
	N::Vector{T} where T<:Integer)

	Cloud(d,N,Array{Float64,2}(undef,Npoints(N),d))
end

function Cscalar(cl::Cloud)
	Cscalar(cl, Vector{ComplexF64}(undef,Npoints(cl)))
end

Npoints(N::Vector{T} where T<:Integer) = prod(N)
Npoints(cl::Cloud)		= Npoints(cl.N) 
Npoints(cs::Cscalar)	= Npoints(cs.cl)
Npoints(scl::SubCloud)	= Npoints(scl.N)
Npoints(scs::SubCscalar)= Npoints(scs.N)

function init!(
	cl::Cloud,
	L::Vector{T} where T<:AbstractFloat)

	if cl.dim==1
		for i=1:Npoints(cl)
			cl.points[i,1] = L[1]*(-1 + 2*(i-1)/cl.N[1])
		end
	elseif cl.dim==2
		for i=1:Npoints(cl)
			cl.points[i,1] = L[1]*(-1 + 2*rem(i-1,cl.N[1])/cl.N[1])
			cl.points[i,2] = L[2]*(-1 + 2*div(i-1,cl.N[1])/cl.N[2])
		end
	elseif cl.dim==3
		for i=1:Npoints(cl)
			cl.points[i,1] = L[1]*(-1 + 2*rem(i-1,cl.N[1])/cl.N[1])
			cl.points[i,2] = L[2]*(-1 + 2*rem(div(i-1,cl.N[1]),cl.N[2])/cl.N[2])
			cl.points[i,3] = L[3]*(-1 + 2*div(i-1,cl.N[1]*cl.N[2])/cl.N[3])
		end
	else
		throw(ArgumentError("cloud dimension must be 1,2, or 3"))
	end
end

function get_originindex(N::Vector{T} where T<:Integer)
	if length(N)==1
		return div(N[1],2)
	else
		return div(prod(N),2) + get_originindex(N[1:end-1])
	end
end
get_originindex(cl::Cloud)	= get_originindex(cl.N)
get_originindex(cs::Cscalar)= get_originindex(cs.cl.N)

# sample a function, f, on a cloud as a complex scalar field
function sample!(cs::Cscalar, f::Function)
	for i=1:Npoints(cs.cl)
		cs.values[i] = f(cs.cl.points[i,1:end])
	end
end

# printers
function tovtk(cl::Cloud, filename::AbstractString)
	f = open(filename, "w")
	writepreamble(f, "Point cloud")
	writepoints(f, cl.dim, Npoints(cl), cl.points)
	writeverts(f,Npoints(cl))
	close(f)
end

function tovtk(scl::SubCloud, filename::AbstractString)
	f = open(filename, "w")
	writepreamble(f, "Point (sub) cloud")
	writepoints(f, scl.cl.dim, Npoints(scl), scl.cl.points[scl.verts,1:end])
	writeverts(f,Npoints(scl))
	close(f)
end

function tovtk(cs::Cscalar, filename::AbstractString)
	f = open(filename, "w")
	writepreamble(f, "Complex scalar field")
	writepoints(f, cs.cl.dim, Npoints(cs.cl), cs.cl.points)
	writeverts(f,Npoints(cs.cl))
	print(f, "CELL_DATA "*string(Npoints(cs.cl))*"\n")
	writescalars(f,"real",Npoints(cs.cl),real.(cs.values))
	writescalars(f,"imag",Npoints(cs.cl),imag.(cs.values))
	close(f)
end

function tovtk(scs::SubCscalar, filename::AbstractString)
	f = open(filename, "w")
	writepreamble(f, "Complex (sub) scalar field")
	writepoints(f, scs.cs.cl.dim, Npoints(scs), scs.cs.cl.points[scs.verts,1:end])
	writeverts(f,Npoints(scs))
	print(f, "CELL_DATA "*string(Npoints(scs))*"\n")
	writescalars(f,"real",Npoints(scs),real.(scs.cs.values[scs.verts]))
	writescalars(f,"imag",Npoints(scs),imag.(scs.cs.values[scs.verts]))
	close(f)
end

function writepreamble(f::IO, label::AbstractString)
	preamble = """
	# vtk DataFile Version 4.0
	"""*label*"""\n
	ASCII
	DATASET UNSTRUCTURED_GRID
	"""
	print(f,preamble)
end

function writepoints(
	f::IO, 
	DIMENSION::Integer, 
	N::Integer, 
	points::Array{Float64,2})

	print(f, "POINTS $N float\n")
	for i=1:N
        for j=1:DIMENSION
            print(f, "$(points[i,j]) ")
        end
        for j=DIMENSION+1:3
            print(f, "0.0 ")
        end
        print(f, "\n")
    end
end

function writeverts(f::IO, N::Integer)
    print(f,"CELLS $N $(2N)\n")
    for i=1:N
        print(f, "1 $(i-1)\n")
    end
    print(f,"CELL_TYPES $N\n")
    for i=1:N
        print(f, "1\n")
    end
end

function writescalars(
        f::IO,
        name::AbstractString,
        N::Integer,
        scalars::Vector{Float64}
    )
    print(f, "SCALARS "*name*" float\nLOOKUP_TABLE default\n")
    for i=1:N
        print(f, "$(scalars[i])\n")
    end
end

function tocsv(scs::SubCscalar, filename)
	f = open(filename, "w")
	if scs.cs.cl.dim==3
		print(f, "x,y,z,Re,Im\n")
		for i=1:Npoints(scs)
			print(f,
				"$(scs.cs.cl.points[scs.verts[i],1]),"*
				"$(scs.cs.cl.points[scs.verts[i],2]),"*
				"$(scs.cs.cl.points[scs.verts[i],3]),"*
				"$(real(scs.cs.values[scs.verts[i]])),"*
				"$(imag(scs.cs.values[scs.verts[i]]))\n")
		end
	elseif scs.cs.cl.dim==2
		print(f, "x,y,Re,Im\n")
		for i=1:Npoints(scs)
			print(f,
				"$(scs.cs.cl.points[scs.verts[i],1]),"*
				"$(scs.cs.cl.points[scs.verts[i],2]),"*
				"$(real(scs.cs.values[scs.verts[i]])),"*
				"$(imag(scs.cs.values[scs.verts[i]]))\n")
		end
	elseif scs.cs.cl.dim==1
		print(f, "x,Re,Im\n")
		for i=1:Npoints(scs)
			print(f,
				"$(scs.cs.cl.points[scs.verts[i],1]),"*
				"$(real(scs.cs.values[scs.verts[i]])),"*
				"$(imag(scs.cs.values[scs.verts[i]]))\n")
		end
	else
		throw(ArgumentError("invalid number of dimensions specified."))
	end
	close(f)
end	