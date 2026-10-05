function eval_val_arrays!(V::Validation, states, polASG, qASG, vDefaultASG, vRepayASG)
    @inbounds Threads.@threads for I in eachindex(states)
        s = states[I]
        V.polVal[I] = polASG(s)
        V.qVal[I] = qASG(s)
        V.vRVal[I] = vRepayASG(s)
    end
    @inbounds for iY in 1:ny
        V.vDVal[iY] = vDefaultASG(SVector(ygrid[iY]))
    end
    return V
end

function create_init_val(states, polASG, qASG, vDefaultASG, vRepayASG)
    V = Validation(zeros(nb, ny),
                   zeros(nb, ny),
                   zeros(ny),
                   zeros(nb, ny))
    eval_val_arrays!(V, states, polASG, qASG, vDefaultASG, vRepayASG)
    return V
end

function update_and_diff!(V::Validation, states, polASG, qASG, vDefaultASG, vRepayASG)
    N = length(states)
    part = collect(Iterators.partition(1:N, cld(N, Threads.nthreads())))
    local_diffs = Vector{NTuple{3,Float64}}(undef, length(part))

    Threads.@threads :static for c in eachindex(part)
        dp = dq = dvR = 0.0   # registers, no false sharing
        @inbounds for I in part[c]
            s = states[I]
            np = polASG(s); 
            nq = qASG(s);   
            nvR = vRepayASG(s)
            dp = max(dp, abs(np - V.polVal[I]))
            dq = max(dq, abs(nq - V.qVal[I]))
            dvR = max(dvR, abs(nvR - V.vRVal[I]))
            V.polVal[I]  = np;
            V.qVal[I]    = nq;
            V.vRVal[I]   = nvR
        end
        local_diffs[c] = (dp, dq, dvR)
    end

    dp = maximum(x->x[1], local_diffs);
    dq = maximum(x->x[2], local_diffs);
    dvR = maximum(x->x[3], local_diffs);

    dvD = 0.0
    @inbounds for iY in 1:ny
        nvD = vDefaultASG(SVector(ygrid[iY]))
        dvD = max(dvD, abs(nvD - V.vDVal[iY]))
        V.vDVal[iY] = nvD
    end
    return dp, dq, dvD, dvR
end
