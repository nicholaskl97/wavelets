# Audio Compression via Discrete Wavelet Transform

This is my final project for the MIT course 18.335, Introduction to Numerical Methods, written in Julia.
The paper (stored in this repository as a PDF) reviews the discrete wavelet transform (DWT) and explains this code.
The code is designed to compute Daubechies wavelet coefficients of arbitrary order, apply a DWT to an audio signal, apply truncated wavelet approximation to compress the signal with minimal losses, and decompress/transform to recover an audio signal.
Furthermore, it compares truncated wavelet approximation (using several sets of Daubechies coefficients) to truncated approximation using discrete cosine series.
The results, shown in the paper, are that in some cases DWT can get better compression for the same mean-squared-error in the decompressed signal than DCT.
While mean-squared-error was used as the error metric, a human listener will also often prefer the DWT-compressed signal, possibly since our ears rely partially on a DWT (see the paper, which gives a reference for this claim).

The test audio used in this project (Am I Blue.wav) is Kevin Conroy's recording of *Am I Blue*, available on his [website](http://www.realkevinconroy.com/).
