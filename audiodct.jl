using FFTW

function compressaudiodct(channels, threshold)
    # Split channels
    left = channels[:,1]
    right = channels[:,2]

    # Perfom transforms
    fleft = dct(left)
    fright = dct(right)

    # Compress
    tleft = threshold * maximum(abs.(fleft))
    ileft = abs.(fleft) .> tleft
    tright = threshold * maximum(abs.(fright))
    iright = abs.(fright) .> tright

    # Generate return values
    fleft = fleft[ileft]
    fright = fright[iright]
    s = sizeof(fleft) + sizeof(fright) + sizeof(ileft) + sizeof(iright) +
        sizeof(length(left))
    return ((fleft, ileft, fright, iright, length(left)), s)
end

function decompressaudiodct(compression)
    # Unpack compression
    cleft, ileft, cright, iright, len = compression

    # Decompress left channel
    fleft = zeros(length(ileft))
    j = 1
    for i = 1:length(fleft)
        if ileft[i]
            fleft[i] = cleft[j]
            j += 1
        end
    end

    # Decompress right channel
    fright = zeros(length(iright))
    j = 1
    for i = 1:length(fright)
        if iright[i]
            fright[i] = cright[j]
            j += 1
        end
    end

    # Perform inverse transforms
    left = idct(fleft)
    right = idct(fright)

    return [left[1:len] right[1:len]]
end
