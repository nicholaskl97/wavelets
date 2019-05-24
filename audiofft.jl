using FFTW

function compressaudiofft(channels, threshold)
    # Split channels
    left = channels[:,1]
    right = channels[:,2]

    # Perfom transforms
    fleft = fft(left)
    fright = fft(right)

    # Compress
    tleft = threshold * maximum(abs.(fleft))
    ileft = abs.(fleft) .> tleft
    tright = threshold * maximum(abs.(fright))
    iright = abs.(fright) .> tright
    return (fleft[ileft], ileft, fright[iright], iright, length(left))
end

function decompressaudiofft(compression)
    # Unpack compression
    cleft, ileft, cright, iright, len = compression

    # Decompress left channel
    fleft = Array{Complex{Float64},1}(undef, length(ileft))
    j = 1
    for i = 1:length(fleft)
        if ileft[i]
            fleft[i] = cleft[j]
            j += 1
        end
    end

    # Decompress right channel
    fright = Array{Complex{Float64},1}(undef, length(iright))
    j = 1
    for i = 1:length(fright)
        if iright[i]
            fright[i] = cright[j]
            j += 1
        end
    end

    # Perform inverse transforms
    left = ifft(fleft)
    right = ifft(fright)

    return real.([left[1:len] right[1:len]])
end
