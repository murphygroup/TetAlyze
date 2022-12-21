# TetAlyze
## Reproduce figures in the paper
1. Setup the paths by running `matlab setup.m`.
2. To reproduce figures in the paper, run the corresponding MATLAB scripts, e.g. 	`matlab fig*.m`.


## Image processing pipeline
0. Raw images used in this project are accessible in [image repo](https://github.com/murphygroup/PearsonGroupImages).
1. Setup the paths by running `matlab setup.m`.
2. Edit the `filename` variable in the `img_proc_pipeline.m` script with the right image path.
3. `rejection_threshold` refers to the distance threshold. The BB with a distance to the fitted
cell outline larger than this threshold will be rejected. The purpose is to reject the noise
puncta inside the cell. The value can be 0.5 (strong rejection) ~ 2 (weak rejection) for WT
cells without big concaves. 
Reduce this value can reject more Noisy puncta inside the cell, which can prevent rows
produced inside the cell.
4. `minRowLength` refer s to the threshold of BB length. The row shorter than this threshold
will be discarded.
5. `minBBsInRow` refers to the threshold of minimum number of BBs in a row. The row with
no greater BBs than this value will be discarded.
6. Consider the importance of the correct row indexing, the code allows the user to correct
the row indexing if necessary. When the Figure 2 window firstly promoted, please enter
the current index of the real first index and press `Enter` to correct the index. If the current
index is correct, enter nothing and hit Enter to continue. 

Image processing was done on a desktop computer with a Core i9-10900K CPU @ 3.70GHz and 64 GB of RAM, running a Windows 11 operating system. The software was tested on MATLAB 2020b.

Required toolbox:
- Bioinformatics Toolbox

- Computer Vision System Toolbox

- Curve Fitting Toolbox

- Image Processing Toolbox

- Signal Processing Toolbox

- Statistics and Machine Learning Toolbox

- Wavelet Toolbox

Note that we include MATLAB Bio-Formats and Locally Weighted Polynomials toolbox packages in our repo here, and we thank the authors for developing those tools.

Images were stored in uint16 TIFF format.

For an image with 500x500x75 pixels, it typically takes about 2 minutes to finish image processing.
