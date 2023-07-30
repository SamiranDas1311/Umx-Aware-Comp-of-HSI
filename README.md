# Umx-Aware-Comp-of-HSI
Unmixing Aware Compression of Hyperspectral Image
First, compute the number of endmembers in hyperspectral data using MixedNormEval.m
Next, perform library pruning by the CovarianceSimMetric.m
Next, compute the abundance by sparse inversion
Find out the compressible abundance maps by H_Mteric.m
Compress these compressible abundance maps by OrthParafac.m
Reconstruct the Hyperspectral Data by Recon.m
