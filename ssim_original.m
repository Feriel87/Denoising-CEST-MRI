function mssim = ssim_original(Org_data,Estim_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ssim_original is a function to calculate the SSIM index beween original and denoised data.
% Function Interface:
%     mssim = ssim_original(Org_data,Estim_data)
% Input Arguments
%     Org_data : Original Data
%     Estim_data : Denoised Data
% Output: 
%     mssim: the calculated ssim value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n] = size(Org_data);
L=255; 
K=[0.01 0.03];
C1 = (K(1)*L)^2;
C2 = (K(2)*L)^2;
ssim_map = zeros(m,n);

IM_original = [Org_data NaN([m 7]); NaN([7 m+7])]; 
IM_Estimate = [Estim_data NaN([m 7]); NaN([7 m+7])];
        for ix=1:m
            for iy =1:n
                
                data_estimate = IM_Estimate(ix:ix+7,iy:iy+7); 
                data_estimate_sq =  data_estimate(:).*data_estimate(:); 
                
                data_org = IM_original(ix:ix+7,iy:iy+7); 
                data_org_sq =  data_org(:).* data_org(:);
                
                data_org_estimate = data_estimate(:).*data_org(:); 
                
                mu1=mean(data_estimate(:),'omitnan');
                mu2=mean(data_org(:),'omitnan');
                
                mu1_sq = mu1.*mu1;
                mu2_sq = mu2.*mu2;
                mu1_mu2 = mu1.*mu2;
                
                sigma1_sq = mean(data_estimate_sq(:),'omitnan')- mu1_sq;
                sigma2_sq = mean(data_org_sq(:),'omitnan')- mu2_sq;
                sigma12 = mean(data_org_estimate(:),'omitnan') - mu1_mu2;

                
                ssim_map(ix,iy) = ((2*mu1_mu2 + C1).*(2*sigma12 + C2))./((mu1_sq + mu2_sq + C1).*(sigma1_sq + sigma2_sq + C2));
            end
        end
        
        mssim = mean2(ssim_map);
end