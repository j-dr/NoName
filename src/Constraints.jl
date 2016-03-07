include("ParticleMesh.jl")

function loadConstraintVector(filename)
    #read in a file containing the constraint
    
    constr = readdlm(filename)
    constr[1:end-1,:], constr[end,:]

end

function idConstraint(ids, dims)

    Npart = prod(dims)
    H = zeros(dims)
    H[ids] = 1/Npart

end

function evalConstr(C, f)

    cf = zeros(f)
    ndims = copy(collect(dims))
    ndims[1] = div(dims[1],2)+1
    
    lsub = dims[1]*dx
    dk = 2*pi/lsub
    d3k = dk*dk*dk

    Hi = conj(rfft(C[i]))
    for ix in 1:ndims[1], iy in 1:dims[2], iz in 1:dims[3]
        
        v = Hi[ix,iy,iz]*f[ix,iy,iz]*d3k
        if (ix>0) & (ix<ndims[1])
            v*=2
        end
        
        cf[ix,iy,iz] += v
    end

    cf
end

    

end

function constrNoise(C, f, Cij, ec, ecp, dims, dx, ps)

    nconst = length(C)
    ndims = copy(collect(dims))
    ndims[1] = div(dims[1],2)+1
    
    lsub = dims[1]*dx
    dk = 2*pi/lsub
    d3k = dk*dk*dk

    for i in 1:nconst
        Hi = conj(rfft(C[i]))
        for j in 1:i
            Hj = rfft(C[j])            
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
        
    f

end


function icovConstr(C, dims, dx, ps)

    nconst = length(C)
    Cij = zeros((nconst, nconst))
    ndims = copy(collect(dims))
    ndims[1] = div(dims[1],2)+1
    
    lsub = dims[1]*dx
    dk = 2*pi/lsub
    d3k = dk*dk*dk

    for i in 1:nconst
        Hi = conj(rfft(C[i]))
        for j in 1:i
            c1 = 0.0
            c2 = 0.0
            Hj = rfft(C[j])
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
    


