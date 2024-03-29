% CEST denoising demo for denoising CEST MRI images. 
% Several denoising filters can be applied to CEST MRI images such as our developed NLmCED filter, 
% BM3D filter, Gaussian filter and Smoothing Cubic Spline method. 
% Different Z-spectra datasets (ground-truth and Rician-noise corrupted) have been provided to test the algorithm such as: 
% 3 synthetic datasets simulated including water, iopamidol and semisolid component (ssMT) pools 
% 2 in-vitro phantoms with Iopamidol Z-spectra at several concentration and pH, acquired at 7T 
% 1 in-vivo data with Z-spectra before (Pre) and after (Post) Iopamidol injection in a murine tumor model.

% **************************************  DATA *************************************************************%

% dataset_1.mat: synthetic dataset #1 (Iopamidol + water pools, chessbox shape)                                          
%             - Data: synthetic noise-free image,                                       
%             - noisyimage1, noisyimage3, noisyimage5, noisyimage7: synthetic images corrupted with 1,3,5 and 7% of Rician noise, respectively.                 
%             - iter1,iter3, iter5, iter7: number of iteration for NLmcED filter to dnoise noisyimage1, noisyimage3, noisyimage5 and noisyimage7 respectively.    
%             - wind: gaussian window size for NLmCED filter (by default =1 for invitro and invivo data)
%             - sig1 ,sig3, sig5, sig7: BM3D parameter to denoise noisyimage1, noisyimage3, noisyimage5 and noisyimage7, respectively.                                                           
%             - reg_spectrum: Smoothing Cubic spline parameter used for denoising noisy images                               
%             - x: frequency offsets                                                                    

% dataset_2.mat: synthetic dataset #2 (Iopamidol + water pools, circular shape)        
%             - Data: synthetic noise-free image,                                       
%             - noisyimage1, noisyimage3, noisyimage5, noisyimage7: synthetic images corrupted with 1,3,5 and 7% of Rician noise respectively.                 
%             - iter1,iter3, iter5, iter7: number of iteration for NLmcED filter to denoise noisyimage1, noisyimage3, noisyimage5 and noisyimage7 respectively.   
%             - wind : gaussian window size for NLmCED filter (by default =1 for invitro and invivo data)
%             - sig1 ,sig3, sig5, sig7: BM3D parameter to denoise noisyimage1, noisyimage3, noisyimage5 and noisyimage7 respectively.                                                           
%             - reg_spectrum: Smoothing Cubic spline parameter used for denoising noisy images                               
%             - x: frequency offsets                  
%             - mask_slice: mask for the full ROI in image                                   
                                                                                                                                                                                                                 
% dataset_3.mat: synthetic dataset #3 (Iopamidol + water + ssMT pools, four tumor-derived shapes, Z-spectra as raw Signal Intensities)  
%             - Data_Pre: synthetic data pre injection data.                
%             - Data_Post: synthetic data post injection data.             
%             - noisypre1, noisypre3, noisypre5, noisypre7: Data_Pre corrupted with 1,3,5 and 7% of Rician noise respectively. 
%             - noisypost1, noisypost3, noisypost5, noisypost7: Data_Post corrupted with 1,3,5 and 7% of Rician noise respectively. 
%             - iter1,iter3,iter5, iter7: number of iteration for NLmcED filter to dnoise corrupted data with 1,3,5 and 7% of Rician noise respectively.  
%             - wind : gaussian window size for NLmCED filter (by default =1 for invitro and invivo data)
%             - sig1 ,sig3, sig5, sig7: BM3D parameter to denoise corrupted data with 1,3,5 and 7% of Rician noise respectively                                                                                                                                 
%             - reg_spectrum: Smoothing Cubic spline parameter used for denoising noisy images(same parameter used for all data).                              
%             - x: frequency offsets   
%             - mask_slice: mask for the full ROI in image                                    
%             - mask_slice_roi: mask for each roi in the full data.  

% invitro_1.mat: in-vitro phantom #1                                                                        
%             - Data: Original noisy data with several pH and Concentration.                        
%             - mask_slice: mask for the full ROI in image                                                    
%             - mask_slice_roi: mask for each roi in the full data.                                              
%             - image: morphological image.                                                         
%             - iter: number of iteration for NLmcED method                                                      
%             - sig: BM3D parameter                                                                              
%             - reg_spectrum: Smoothing Cubic spline parameter (same parameter used for all data).                                                  
%             - x: offsets    

% invitro_2.mat : in-vitro phantom #2                                                                        
%             - Data: Original noisy data with several Concentration at pH = 6.7.
%             - mask_slice: mask for the full ROI in image                                    
%             - mask_slice_roi: mask for each roi in the full data.                           
%             - image: morphological image, size: 64x64                                       
%             - iter: number of iteration for NLmcED method                                   
%             - sig: BM3D parameter                                                           
%             - reg_spectrum: Smoothing Cubic spline parameter (same parameter used for all data).                                
%             - x: frequency offsets   

% invivo_1.mat : in-vivo data                                                                   
%             - Data_Pre: In-vivo pre injection data                          
%             - Data_Post: In-vivo post injection data                      
%             - mask_slice: mask for the full ROI in image                                    
%             - mask_slice_roi: mask for each roi in the full data                            
%             - image: morphological image                                                    
%             - iter: number of iteration for NLmcED filter to denoise  Data_Pre and Data_Post respectively.                                
%             - sig: BM3D parameters to denoise  Data_Pre and Data_Post respectively.                                                           
%             - reg_spectrum: Smoothing Cubic spline parameter(same parameter used for all data).                               
%             - x: frequency offsets 


 
%% Example denoising 1% of rician noise 
load('dataset_1.mat')
[m,n,p]    = size(Data); % Size of the original data

%****** Applying NLmCED filter for denoising 1% of rician noise ******%
%%% rho and alpha are the optimised parameters for the CEST data (synthetic, invitro and invivo); 
%%%  If you are using another modality you can chage these parameter to optimise the denoising process 
rho = 0.001;
alpha = 0.001;
denoisedNLmCED = NLmCED(noisyimage1,iter1,rho, alpha,wind);

%****** Applying BM3D filter for denoising 1% of rician noise ******%
denoisedBM3D = zeros(m,n,p);
for i =1:p
    [~, out] = BM3D(1, noisyimage1(:,:,i),sig1);
    denoisedBM3D(:,:,i)= out;
end

%****** Applying Gaussian filter for denoising 1% of rician noise ******%
denoisedGaussian = zeros(m,n,p);
for i =1:p
    denoisedGaussian(:,:,i) = imgaussfilt(noisyimage1(:,:,i),'FilterSize',7);
end

%****** Applying Smoothing Cubic Spline for denoising 1% of rician noise ******%
denoisedCubicSS = zeros(m,n,p);
for i=1:m
    for j=1:n
        noisy=squeeze(noisyimage1(i,j,:)); % z-spectrum
        cs = csaps(x,noisy,reg_spectrum,[]);
        [minval_mean,min_absc_mean]=fnmin(cs);
        denoisedCubicSS(i,j,:) = fnval(cs,x+min_absc_mean);
        clear minval_mean min_absc_mean cs noisy
    end
end


%% ******************************* Example of calculation of CEST Contrast and pH ******************************** %%
% Example calculation of CEST Contrast and pH for free noise data 
% with only Water+Iopamidol (dataset_1.mat) and for the denoised dataset with NLmCED method
% %Original data
[ST1_Org, ST2_Org] = contrastCEST(Data,x); 
pH_Org = pH_SyntheticDataset(ST1_Org, ST2_Org); 
% %Denoised data
[ST1_denoised, ST2_denoised] = contrastCEST(denoisedNLmCED,x); 
pH_denoised = pH_SyntheticDataset(ST1_denoised, ST2_denoised); 

% calculate the PSNR index for ST1 (Contrast CEST at 4.2 ppm)
PSNR_ST1 = psnr_original(ST1_Org,ST1_denoised); 
SSIM_ST1 = ssim_original(ST1_Org,ST1_denoised); 

%% *************************************************************************************************************** %%
% In order to test data with two sets (pre and post injection), please follow the example below for testing dataset3. 
% % Example Calculation of DeltaST and pH for free noise data with Water+MT+Iopamidol (dataset_3.mat): 
load('dataset_3.mat')
% NLmCED Filter
rho = 0.001;
alpha = 0.001;
denoisedNLmCED_pre = NLmCED(noisypre1,iter1,rho,alpha,wind);
denoisedNLmCED_post = NLmCED(noisypost1,iter1,rho,alpha,wind);
% DeltaST and pH calculation for Original data
[ST1_Org_Pre, ST2_Org_Pre] = contrastCEST(Data_Pre,x);  
[ST1_Org_Post, ST2_Org_Post] = contrastCEST(Data_Post,x); 
deltaST1_Org = ST1_Org_Post - ST1_Org_Pre; 
deltaST2_Org = ST2_Org_Post - ST2_Org_Pre;
pH_Org = pH_SyntheticDataset(deltaST1_Org, deltaST2_Org); 
% DeltaST and pH calculation for denoised NLmCED
[ST1_den_Pre, ST2_den_Pre] = contrastCEST(denoisedNLmCED_pre,x);  
[ST1_den_Post, ST2_den_Post] = contrastCEST(denoisedNLmCED_post,x); 
deltaST1_nlmced = ST1_den_Post - ST1_den_Pre; 
deltaST2_nlmced = ST2_den_Post - ST2_den_Pre; 
pH_nlmced = pH_SyntheticDataset(deltaST1_nlmced, deltaST2_nlmced); 
% PSNR and SSIM calculation
DeltaST1_psnrval = psnr_original(deltaST1_Org,deltaST1_nlmced); 
DeltaST1_ssImval = ssim_original(deltaST1_Org,deltaST1_nlmced); 

%% ****************************************************************************%
% This an Example to calculate the DeltaST and pH for invivo data (invivo_1.mat) 
% Original data
% load('invivo_1.mat');
[ST1_Org_Pre, ST2_Org_Pre] = contrastCEST(Data_Pre,x);  
[ST1_Org_Post, ST2_Org_Post] = contrastCEST(Data_Post,x); 
deltaST1_Org = ST1_Org_Post - ST1_Org_Pre; 
deltaST2_Org = ST2_Org_Post - ST2_Org_Pre;
deltaST_Ratio = deltaST1_Org./deltaST2_Org;
pH_Org = pH_InVivo(deltaST1_Org, deltaST2_Org); 
% Denoised data using NLmCED method
rho = 0.001;
alpha = 0.001;
wind =2; 
denoisedNLmCED_pre = NLmCED(Data_Pre,iter,rho, alpha,wind);
denoisedNLmCED_post = NLmCED(Data_Post,iter,rho,alpha,wind);
% DeltaST and pH calculation for denoised NLmCED
[ST1_den_Pre, ST2_den_Pre] = contrastCEST(denoisedNLmCED_pre,x);  
[ST1_den_Post, ST2_den_Post] = contrastCEST(denoisedNLmCED_post,x); 
deltaST1_nlmced = ST1_den_Post - ST1_den_Pre; 
deltaST2_nlmced = ST2_den_Post - ST2_den_Pre; 
pH_nlmced = pH_InVivo(deltaST1_nlmced, deltaST2_nlmced); 

