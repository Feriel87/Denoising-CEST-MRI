function data_pH =  pH_SyntheticDataset(ST1_data, ST2_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pH_SyntheticDataset is a function to calculate the pH values in the synthetic datasets.
% Function Interface:
%     data_pH =  pH_SyntheticDataset(ST1_data, ST2_data)
% Input Arguments
%     ST1_data : Contrast CEST at 4.2 ppm of iopamidol 
%     ST2_data : Contrast CEST at 5.5 ppm of iopamidol 
% Output: 
%     data_pH: the calculated pH values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n] = size(ST1_data);
data_pH = zeros(m,n); 
 p1 = -6.817; p2 = 15.76; p3 = -12.64; p4 = 10.45;
        DeltaDiv_data = ST2_data./ ST1_data;
        DeltaDiv_data(isnan(DeltaDiv_data))= 0;
        for i=1:m
            for j=1:n
                % Free noise data (ground truth data)
                if (ST1_data(i,j) > 1 && ST2_data(i,j) > 1 )
                    data_pH(i,j)=(p1.*DeltaDiv_data(i,j).^3)+(p2.*DeltaDiv_data(i,j).^2)+(p3.*DeltaDiv_data(i,j))+ p4;
                else
                    data_pH(i,j)=0;
                end
            end
        end
end