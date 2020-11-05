clear all
clc

org=imread('Pelvix_Orig.tif');
nous=imread('Pelvix_Nous.tif'); 
NLM=imread('Pelvix_NLM_01.tif'); 
CED=imread('Pelvix_CED.tif'); 
TV=imread('Pelvix_TV.tif'); 
 
diff_nous = imabsdiff(org,nous);
diff_NLM = imabsdiff(org,NLM);
diff_CED = imabsdiff(org,CED);
diff_TV = imabsdiff(org,TV);

figure(1)
imshow(diff_nous)
figure(2)
imshow(diff_NLM)
figure(3)
imshow(diff_CED)
figure(4)
imshow(diff_TV)

%% 

% [D, META] =nrrdread('5157-005-01_10_T1.nrrd');
%  viewer3d(D)