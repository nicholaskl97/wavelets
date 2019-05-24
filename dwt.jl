using NLsolve

# Precondition: ∃ J ∈ ℤ: length(f) = 2^J
# Returns: ((a0, b0, b1, ...), length(a0))
function mallatdecomp(f, L, H, p)
    out = []
    a = f
    while length(a) >= 2p
        out = [H(a); out]
        a = L(a)
    end
    out = [a; out]
    return (out, length(a))
end

# Precondition: ∃ J ∈ ℤ: length(f) = 2^J,
#   length(a0) = length(b0),
#   length(b_{i+1}) = 2 length(b_i)
# Input: Wf = (a0, b0, b1, ...), n = length(a0)
function mallatrecon(Wf, n, iL, iH)
    a = Wf[1:n]
    i = n
    while i + n <= length(Wf)
        b = Wf[(i+1):(i+n)]
        a = iL(a) + iH(b)
        i += n
        n = 2n
    end
    return a
end

# Returns: out[i+1] = c_i for i ∈ {0, 2p-1} for DAUB(2p) wavelets
function daubgen(p)
    function f!(F,c)
        for m = 1:(p-1)
            k = collect((2m):(2p-1))
            F[m] = sum(c[k .+ 1] .* c[k .- 2m .+ 1])
        end
        F[p] = sum(c.^2) - 2
        k = collect(0:(2p - 1))
        for m = 1:p
            F[p+m] = sum((-1).^k .* k.^(m-1) .* c)
        end
    end

    function j!(J, c)
        for m = 1:(p-1)
            J[m, (2p-2m+1):(2p)] = zeros(1, 2m)
            for i = 0:(2p - 2m - 1)
                J[m,i+1] = c[i+1+2m]
            end
            for i = (2m):(2p-1)
                J[m,i+1] += c[i+1-2m]
            end
        end
        J[p,:] = 2 * c
        k = collect(0:(2p - 1))
        for m = 1:p
            J[p+m, :] = (-1).^k .* k.^(m-1)
        end
    end

    r = nlsolve(f!, j!, ones(2p))
    if converged(r)
        return r.zero
    end
    throw(Exception)
end

function filtergen(p)
    c = daubgen(p)

    function L(f)
        n = length(f)
        out = zeros(n >> 1)
        for i = 1:(n >> 1)
            s = 0
            for j = (2(i-p) + 1):(2i)
                s += f[(n+j-1)%n+1] * c[(2i - j + 1)]
            end
            out[i] = s / 2
        end
        return out
    end

    function iL(g)
        nhalf = length(g)
        n = 2 * nhalf
        out = zeros(n)
        for j = 1:n
            s = 0
            for i = ((j+1)>>1):(p-1 + (j+1)>>1)
                s += g[(i-1)%nhalf+1] * c[(2i - j + 1)]
            end
            out[j] = s
        end
        return out
    end

    function H(f)
        n = length(f)
        out = zeros(n >> 1)
        for i = 1:(n >> 1)
            s = 0
            for j = (2i-1):(2(p+i-1))
                s += f[(j-1)%n+1] * (-1)^(j+1) * c[(j -2i + 2)]
            end
            out[i] = s / 2
        end
        return out
    end

    function iH(g)
        nhalf = length(g)
        n = 2 * nhalf
        out = zeros(n)
        for j = 1:n
            s = 0
            for i = ((j+1)>>1 + 1 - p):((j+1)>>1)
                s += g[(nhalf+i-1)%nhalf+1] * (-1)^(j+1) * c[(j - 2i + 2)]
            end
            out[j] = s
        end
        return out
    end

    return (L, H, iL, iH)
end
