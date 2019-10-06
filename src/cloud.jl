#
#	Types
#
mutable struct Cloud
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

function toint(A::Vector{Int}, N::Vector{Int})
	if length(N)==1
		return A[1]+1
	else
		return A[end]*prod(N[1:end-1]) + toint(A[1:end-1],N[1:end-1])
	end
end

function tovec(i::Int, N::Vector{Int})
	if length(N)==1
		return [rem(i-1,N[1])]
	else
		return vcat(tovec(i,N[1:end-1]),[rem(div(i-1,prod(N[1:end-1])),N[end])])
	end
end

function getverts(cl::Cloud, axis::Integer)

	axis <= dim(cl) || throw(ArgumentError("your cloud is $(dim(cl))-dimensional but you are asking to slice along axis $(axis)"))

	A = div.(cl.N,2)
	A[axis] = 0
	verts = [toint(A,cl.N)]

	for i=1:cl.N[axis]-1
		A[axis] += 1
		verts = vcat(verts,toint(A,cl.N))
	end

	N = ones(Int, dim(cl))
	N[axis] = cl.N[axis]

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

	Cloud(N,Array{Float64,2}(undef,Npoints(N),d))
end

function Cscalar(cl::Cloud)
	Cscalar(cl, Vector{ComplexF64}(undef,Npoints(cl)))
end

dim(N::Vector{T} where T<: Integer) = length(N)
dim(cl::Cloud)   	= dim(cl.N)
dim(scl::SubCloud) 	= dim(scl.cl.N)
dim(cs::Cscalar) 	= dim(cs.cl.N)
dim(scs::SubCscalar)= dim(scs.cs.cl.N)

Npoints(N::Vector{T} where T<:Integer) = prod(N)
Npoints(cl::Cloud)		= Npoints(cl.N) 
Npoints(cs::Cscalar)	= Npoints(cs.cl)
Npoints(scl::SubCloud)	= Npoints(scl.N)
Npoints(scs::SubCscalar)= Npoints(scs.N)

function init!(
	cl::Cloud,
	L::Vector{T} where T<:AbstractFloat)

	for i=1:Npoints(cl)
		A = tovec(i, cl.N)
		@. cl.points[i,1:end] = L*(-1 + 2*A/cl.N)
	end
end

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
	writepoints(f, dim(cl), Npoints(cl), cl.points)
	writeverts(f,Npoints(cl))
	close(f)
end

function tovtk(scl::SubCloud, filename::AbstractString)
	f = open(filename, "w")
	writepreamble(f, "Point (sub) cloud")
	writepoints(f, dim(scl), Npoints(scl), scl.cl.points[scl.verts,1:end])
	writeverts(f,Npoints(scl))
	close(f)
end

function tovtk(cs::Cscalar, filename::AbstractString)
	f = open(filename, "w")
	writepreamble(f, "Complex scalar field")
	writepoints(f, dim(cs), Npoints(cs.cl), cs.cl.points)
	writeverts(f,Npoints(cs.cl))
	print(f, "CELL_DATA "*string(Npoints(cs.cl))*"\n")
	writescalars(f,"real",Npoints(cs.cl),real.(cs.values))
	writescalars(f,"imag",Npoints(cs.cl),imag.(cs.values))
	close(f)
end

function tovtk(scs::SubCscalar, filename::AbstractString)
	f = open(filename, "w")
	writepreamble(f, "Complex (sub) scalar field")
	writepoints(f, dim(scs), Npoints(scs), scs.cs.cl.points[scs.verts,1:end])
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
	if dim(scs)==3
		print(f, "x,y,z,Re,Im\n")
		for i=1:Npoints(scs)
			print(f,
				"$(scs.cs.cl.points[scs.verts[i],1]),"*
				"$(scs.cs.cl.points[scs.verts[i],2]),"*
				"$(scs.cs.cl.points[scs.verts[i],3]),"*
				"$(real(scs.cs.values[scs.verts[i]])),"*
				"$(imag(scs.cs.values[scs.verts[i]]))\n")
		end
	elseif dim(scs)==2
		print(f, "x,y,Re,Im\n")
		for i=1:Npoints(scs)
			print(f,
				"$(scs.cs.cl.points[scs.verts[i],1]),"*
				"$(scs.cs.cl.points[scs.verts[i],2]),"*
				"$(real(scs.cs.values[scs.verts[i]])),"*
				"$(imag(scs.cs.values[scs.verts[i]]))\n")
		end
	elseif dim(scs)==1
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