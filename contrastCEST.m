function [ST1_data, ST2_data] = contrastCEST(data,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% contrastCEST is a function to calculate the CEST contrast for the in-vivo data for the two pools of Iopamidol (4.2 and 5.5 ppm).
% Function Interface:
%    contrastCEST(data,x)
% Input Arguments
%     data : the normalized data
%      x   : the offsets axis 
% Outputs: 
%     ST1_data : CEST Contrast at 4.2 ppm of iopamidol 
%     ST2_data : CEST Contrast at 5.5 ppm of iopamidol 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n,p] = size(data);
ppm1=4.2; %1st pool of Iopamidol
ppm2=5.5; %2nd pool of Iopamidol
for i=1:m
    for j=1:n
        % Calculate the ST for the free noise data (ground truth data)
        Zspec_data=squeeze(data(i,j,:));
        ST1_data(i,j)=(Zspec_data(find(abs(x-(-ppm1)) < 0.001))-Zspec_data(find(abs(x-(ppm1)) < 0.001)))/Zspec_data(find(abs(x-(-ppm1)) < 0.001))*100;
        ST2_data(i,j)=(Zspec_data(find(abs(x-(-ppm2)) < 0.001))-Zspec_data(find(abs(x-(ppm2)) < 0.001)))/Zspec_data(find(abs(x-(-ppm2)) < 0.001))*100;
    end
end
end
