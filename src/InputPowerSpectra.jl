
type scaleFreePS
    n
    norm
end

type multiplePeaksPS 
    k
    lnδk 
    ps::scaleFreePS # could generalize this further to more PSs
end

type cambPS
    k
    p
end

inputPowerSpectrum = Union{scaleFreePS, multiplePeaksPS, cambPS}

multiplePeaksPS(n,norm,k,lnδk) = multiplePeaksPS(k, lnδk,
                                                 scaleFreePS(n,norm))

function cambPS(filename)
    
    nums = readdlm(filename)
    k = nums[:,1]
    p = nums[:,2]

    cambPS(k, p)

end

function multiplePeaksPS{T<:Real}(n, norm, k::Array, lnδk::T)
    ndelta = zeros(k)
    ndelta[:] = lnδk
    multiplePeaksPS(n, norm, k, ndelta)
end

P(k::Real, PS::scaleFreePS) = PS.norm*convert(Float64,k)^PS.n

function P{T<:Real}(ki::T, PS::cambPS)
    idx = searchsorted(PS.k, ki)
    PS.p[idx.start]
end

function P{T<:Real}(ki::T, PS::multiplePeaksPS)
    ret = 0.
    for i in 1:length(PS.k)
        if abs(ki-PS.k[i])/ki < PS.lnδk[i]/2
            ret = P(ki, PS.ps)
        end
    end
    ret
end

# makes sure we can call all of them also with arrays of k vectors
function P{T<:Array}(k::T, PS::inputPowerSpectrum)
    res = zeros(Float64,size(k))
    for i in 1:length(k)
        res[i] = P(k[i], PS)
    end
    res
end