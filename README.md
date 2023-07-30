# Umx-Aware-Comp-of-HSI
This work performs compression of hyperspectral data using library based unmixing and PARAFAC tensor decomposition
The code requires tensor toolbox version 2.6.
The tensor toolbox is available at: %https://www.sandia.gov/~tgkolda/TensorToolbox/index-2.6.html
Steps
1. Unmixing Aware Compression of Hyperspectral Image
2. First, compute the number of endmembers in hyperspectral data using MixedNormEval.m
3. Next, perform library pruning and sparse inversion by LibUnmJBLD.m
4. Find out the compressible abundance maps by H_Mteric.m
5. Compress these compressible abundance maps by OrthParafac.m
7. Reconstruct the Hyperspectral Data by Recon.m
