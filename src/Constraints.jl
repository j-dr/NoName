#include("ParticleMesh.jl")

function fftfreq(iijk, dims)
    ijk = iijk
    s = div(dims,2)
    k = 2pi .* (ijk-1 - ((ijk .> (s+1)).*dims)) ./dims
end

function loadConstraintVector(filename)
    #read in a file containing the constraint
    
    cfp = open(filename, "r")
    
    #read header containing number of constraints
    #and length of constraint support for each
    hdr = read(cfp, Int64, 2)
    
    #Read value of constraints
    ci = read(cfp, Float64, hdr[1])
    
    #read particle ids to apply constraints to
    pid = read(cfp, Int64, prod(hdr))

    hdr, ci, pid
end

function idConstraint(ids, dims)

    Npart = prod(dims)
    H = zeros(dims)
    H[ids] = 1/Npart

end

function getConstraints(nconstr, csupp, pids)

    dims = conf["ParticleDimensions"]
    Npart = prod(dims)
    alpha = zeros(Npart, nconstr)
    count = 1
    for i in 1:nconstr
        alpha[pids[count:count+csupp[i]-1],i] = 1/Npart
        count += csupp[i]
    end
    alpha

end

function evalConstr(C, f, dx)

    dims = conf["ParticleDimensions"]
    cr = zeros(Complex{Float64}, nconstr)
    ndims = copy(collect(dims))
    ndims[1] = div(dims[1],2)+1
    
    lsub = dims[1]*dx
    dk = 2*pi/lsub
    d3k = dk*dk*dk
    println(dx)
    println(d3k)
    println(lsub)
    println(dims)
    for i in 1:nconstr
        Hi = conj(rfft(reshape(C[:,i], tuple(dims...))))
        for ix in 1:ndims[1], iy in 1:dims[2], iz in 1:dims[3]
            v = Hi[ix,iy,iz]*f[ix,iy,iz]*d3k
            if (ix>0) & (ix<ndims[1])
                v*=2
            end
            cr[i] += v
        end
    end
    println(cr[1].im)
    println(cr[1].re)

    cr = [c.re for c in cr]
end

function constrNoise(C, f, Cij, ec, ecp, dims, dx, ps)

    ndims = copy(collect(dims))
    ndims[1] = div(dims[1],2)+1
    
    lsub = dims[1]*dx
    dk = 2*pi/lsub
    d3k = dk*dk*dk

    for i in 1:nconstr
        Hi = conj(rfft(reshape(C[:,i], tuple(dims...))))
        for j in 1:i
            Hj = rfft(reshape(C[:,j], tuple(dims...)))            
            for ix in 1:ndims[1], iy in 1:dims[2], iz in 1:dims[3]
                kf = fftfreq([ix,iy,iz], dims)
                k = sqrt(dot(kf,kf))
                if k==0
                    continue
                end
                v = Hi[ix,iy,iz]*Cij[i,j]*(ec[i]-ecp[j])*P(k,ps)*d3k
                v += Hj[ix,iy,iz]*Cij[j,i]*(ec[j]-ecp[i])*P(k,ps)*d3k

                if (ix>0) & (ix<ndims[1])
                    v*=2
                end

                f[ix,iy,iz] += v
            end
        end
    end
end


function icovConstr(C, dims, dx, ps)

    Cij = zeros((nconstr, nconstr))
    ndims = copy(collect(dims))
    ndims[1] = div(dims[1],2)+1
    
    lsub = dims[1]*dx
    dk = 2*pi/lsub
    d3k = dk*dk*dk

    for i in 1:nconstr
        Hi = conj(rfft(reshape(C[:,i], tuple(dims...))))
        for j in 1:i
            println(string("i: ",i))
            println(string("j: ",j))
            c1 = 0.0
            c2 = 0.0
            Hj = rfft(reshape(C[:,j], tuple(dims...)))
            for ix in 1:ndims[1], iy in 1:dims[2], iz in 1:dims[3]
                kf = fftfreq([ix,iy,iz], dims)
                k = sqrt(dot(kf,kf))
                if k==0
                    continue
                end
                v = Hi[ix,iy,iz]*Hj[ix,iy,iz]*P(k,ps)*d3k
                if (ix>0) & (ix<ndims[1])
                    v*=2
                end

                c1 += v.re
                c2 += conj(v).re
            end
            Cij[i,j] = c1
            Cij[j,i] = c2
        end
    end
    
    inv(Cij)

end
    
function applyConstraints(c, dx, dims)
    
    hdr, ci, pid = loadConstraintVector(constrfile)
    global nconstr = hdr[1]
    alpha = getConstraints(hdr[1], hdr[2], pid)
    cr = evalConstr(alpha,c,dx)
    Cij = icovConstr(alpha, dims, dx, ps)
    constrNoise(alpha, c, Cij, cr, ci, dims, dx, ps)

end