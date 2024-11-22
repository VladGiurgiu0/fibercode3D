# fibercode3D
Matlab codes to process 3D MART objects containing tracers and fibers.

These codes have been developed in collaboration with Marco de Paoli and Mobin Alipour.

Example input data can be found at: https://owncloud.tuwien.ac.at/index.php/s/pA4J51aZE1WBCQ3
For an example of use of this code see: https://doi.org/10.48550/arXiv.2406.12462

Procedure:
1. Process images with MART in Davis with the option "Store one file for each z-plane" (Reconstruction -> Enable advanced settings).
2. Use "Vlad_loops.m" to set your processing parameters and indicate the folders where your MART is (note: the "p.load" parameter should point to a folder like e.g. "J:\Fibres_spinning_close_wall_Day3\LOOP_Re720_AfterLaserModification\ImgPreproc\TomographicPIV\loop=0\Data\") and where the fiber data should be saved. This code calls "Vlad_main_for_loops.m" which does 4 things: a) discriminates the fibers, b) tracks the fibers, c) models the fibers as polynomials, and d) computes fiber quantities such as spinning and tumbling rates. Note: "loop" means a subfolder containing one data set out of multiple statistically independent sets ("loops") of time-resolved images.
3. (Optional) Use "Vlad_plot_statistics_v5.m" to collect all the data from the loops, to remove outliers (based on proximity to edges of the volume or other quantities), and compute statistics of the fiber quantities.

Notes:
1. "synthethic_fiber_generator" codes can be used to generate 3D fibers rotating at given spinning and tumbling rates to test the main fiber tracking code.


Things to do:
- Renormalize the red vector to have the length.
- Integrate the method to fit the rotation rates like Voth.
- Remove the time-skipping?
- Clean Eigenvector_Tensor. Save principal axis only under one name.
- Opt: Make an option so that when loading MART files binarize them before giving it to be discriminated -> saves RAM and so we can run on multiple CPUs at the same time.
