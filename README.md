### What is this repository for? ###

  Demo_Denoising is a collection of Matlab-based scripts for denoising CEST MRI images. 
  Several filters have been provided for denoising MRI-CEST images, including our developed NLmCED filter and several established filters: BM3D, Gaussian and Smoothing Cubic Spline method.
  
  Several Z-spectra datasets have been provided to test the scripts, including:
  
  3 synthetic datasets by simulating Iopamidol at several concentrations, pH values and amplitudes of the semisolid component including ground-truth data and noisy data after    applying rician noise at several noise levels 
  
  2 in-vitro phantom with Iopamidol at several pH values and concentrations
  
  1 in-vivo datasets with Z-spectra acquired before (PRE) and after (POST) Iopamidol injection in a tumor murine model. 

Copyright (c) 2020-2030   University of Al Manar Tunisia, (ENIT), Tunisia and University of Turin (Unito), Italy

All rights reserved.
This work should only be used for nonprofit purposes.

Authors: Feriel Romdhane

### Installation ###
Unzip CEST_denoising.zip (contains codes and test data) in a folder that is in the MATLAB path.
Demo_CEST_Denoising.m: Demo on how to use the denoising method and data.
contrastCEST.m:        Function to calculate the CEST Contrast for iopamidol.
pH_SyntheticDataset.m: Function to calculate the pH for the synthetic dataset.
pH_InVivo.m:           Function to calculate the pH for the In Vivo data.
psnr_original.m:       Function to calculate the PSNR index.
ssim_original.m:       Function to calculate the SSIM index.
BM3D Filter :          Folder with BM3D package.
NLmCED_Filter:         Folder with NLmCED package.
dataset_1:             Synthetic dataset #1.
invitro_1:             In-vitro phantom #1.
invivo_1:              In-vivo data #1.


### Instructions to download additional data  ###
To test the additional synthetic dataset #3 and the in-vitro phantom #2:

1) Access to the XNAT platform of the Molecular Imaging Center of the University of Torino: http://cim-xnat.unito.it with the following credentials:
username: CESTMRI
password: denoising

2) Click on the project named CEST-MRI-denoising.

3) In the Actions Menu on the right panel, click on Manage Files.

4) A window with a File Manager console will open up, browse in the resource folder named matlab_data to download the data. 



###  Requirements ###

1) MS Windows (32 or 64 bit), Linux (32 bit or 64 bit)
   or Mac OS X (32 or 64 bit)
2) Matlab R2015a or earlier.

###  References ###

[1] F. Romdhane, D. Villano, P. Irrera, L. Consolino, D. L. Longo,
"Improving contrast quantification of MRI-CEST images by applying a denoising approach 
based on a new combination between non-local means filter and anisotropic diffusion tensor,"
in the 14th European Molecular Imaging Meeting (EMIM), Glasgow, UK, March 2019.

[2] K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, "Image 
denoising by sparse 3D transform-domain collaborative filtering," 
IEEE Trans. Image Process., vol. 16, no. 8, August 2007.

[3] DL.Longo, W. Dastru W, G. Digilio, J. Keupp, S. Langereis, 
S. Lanzardo, S. Prestigio, O. Steinbach, E.Terreno, F. Uggeri, S. Aime, 
"Iopamidol as a responsive MRI-chemical exchange saturation transfer 
contrast agent for pH mapping of kidneys: In vivo studies in mice
at 7 T," Magn Reson Med 2011;65(1):202-211.s

### Citation ###
Cite the code: 
[![DOI](https://zenodo.org/badge/310314628.svg)](https://zenodo.org/badge/latestdoi/310314628)

Cite the paper: 
https://doi.org/10.1002/mrm.28676
### Who do I talk to? ###
If you have any comment, suggestion, or question, please do contact Feriel Romdhane at  ferielromdhane@yahoo.fr
