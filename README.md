# fibercode3D
Matlab codes to process 3D MART objects containing tracers and fibers.

These codes have been developed in collaboration with Marco de Paoli and Mobin Alipour.

For an example of use of this code see: https://doi.org/10.48550/arXiv.2406.12462

Procedure:
1. Process images with MART in Davis with the option "Store one file for each z-plane" (Reconstruction -> Enable advanced settings).
2. Use "Vlad_loops.m" to set your processing parameters and indicated the folders where your MART is and where the fiber data should be saved. This code calls "Vlad_main_for_loops.m" which does 4 things: a) discriminates the fibers, b) tracks the fibers, c) models the fibers as polynomials, and d) computes fiber quantities such as spinning and tumbling rates. Note: "loop" means a subfolder containing one data set out of multiple statistically independent sets ("loops") of time-resolved images.
3. (Optional) Use "Vlad_plot_statistics_v5.m" to collect all the data from the loops, to remove outliers (based on proximity to edges of the volume or other quantities), and compute statistics of the fiber quantities.

Notes:
1. "synthethic_fiber_generator" codes can be used to generate 3D fibers rotating at given spinning and tumbling rates to test the main fiber tracking code.


Things to do:
1. Renormalize the red vector to have the length 1.
2. Remove all parameters not needed.
3. Integrate the method to fit the rotation rates like Voth.
4. Create an "example_data" folder.
5. Remove the technique to fit the fiber in the "Plane fitting" way.
6. Remove the std technique for the principal axis.
7. Remove the time-skipping?
8. Clean Eigenvector_Tensor. Save principal axis only under one name.
10. Opt: Make an option so that when loading MART files binarize them before giving it to be discriminated -> saves RAM and so we can run on multiple CPUs at the same time.
11. Opt: Make an option to save keep or delete the fiber images after the code does not need them anymore, so less data is saved.
