include("dwt.jl")

function compressaudio(channels, p, l, h, threshold)
    # Split channels
    left = channels[:,1]
    right = channels[:,2]

    # Pad signals with zeros
    len = length(left)
    J = Integer(ceil(log2(length(left))))
    N = 2^J
    left = [left; zeros(N - len)]
    right = [right; zeros(N - len)]

    # Perfom transforms
    wleft, nstop = mallatdecomp(left,l,h,p)
    wright, _ = mallatdecomp(right,l,h,p)

    # Compress
    tleft = threshold * maximum(abs.(wleft))
    ileft = abs.(wleft) .> tleft
    tright = threshold * maximum(abs.(wright))
    iright = abs.(wright) .> tright

    # Generate return values
    wleft = wleft[ileft]
    wright = wright[iright]
    s = sizeof(wleft) + sizeof(wright) + sizeof(ileft) + sizeof(iright) +
        sizeof(nstop) + sizeof(len)
    return ((wleft, ileft, wright, iright, nstop, len), s)
end

function decompressaudio(compression, p, il, ih)
    # Unpack compression
    cleft, ileft, cright, iright, nstop, len = compression

    # Decompress left channel
    wleft = zeros(length(ileft))
    j = 1
    for i = 1:length(wleft)
        if ileft[i]
            wleft[i] = cleft[j]
            j += 1
        end
    end

    # Decompress right channel
    wright = zeros(length(iright))
    j = 1
    for i = 1:length(wright)
        if iright[i]
            wright[i] = cright[j]
            j += 1
        end
    end

    # Perform inverse transforms
    left = mallatrecon(wleft,nstop,il,ih)
    right = mallatrecon(wright,nstop,il,ih)

    # Unpad
    left = left[1:len]
    right = right[1:len]
    return [left right]
end
