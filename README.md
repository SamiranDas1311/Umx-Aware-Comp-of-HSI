# Umx-Aware-Comp-of-HSI
1. Unmixing Aware Compression of Hyperspectral Image
2. First, compute the number of endmembers in hyperspectral data using MixedNormEval.m
3. Next, perform library pruning by the LibUnmJBLD.m
4. Next, compute the abundance by sparse inversion
5. Find out the compressible abundance maps by H_Mteric.m
6. Compress these compressible abundance maps by OrthParafac.m
7. Reconstruct the Hyperspectral Data by Recon.m
