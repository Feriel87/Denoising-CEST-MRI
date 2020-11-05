function data_pH =  pH_InVivo(ST1_data, ST2_data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pH_InVivo is a function to calculate the pH values for the in-vivo data.
% Function Interface:
%     data_pH =  pH_InVivo(ST1_data, ST2_data)
% Input Arguments
%     ST1_data : CEST Contrast  at 4.2 ppm of iopamidol 
%     ST2_data : CEST Contrast  at 5.5 ppm of iopamidol 
% Output: 
%     data_pH: the calculated pH values 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n] = size(ST1_data);
data_pH = zeros(m,n); 
 p1 = 0.843;  p2 = -3.682; p3 = 5.839; p4 = 3.549;
        DeltaDiv_data = ST1_data./ ST2_data;
        DeltaDiv_data(isnan(DeltaDiv_data))= 0;
        for i=1:m
            for j=1:n
                % Free noise data (ground truth data)
                if (ST1_data(i,j) > 2 && ST2_data(i,j) > 1 )
                    data_pH(i,j)=(p1.*DeltaDiv_data(i,j).^3)+(p2.*DeltaDiv_data(i,j).^2)+(p3.*DeltaDiv_data(i,j))+ p4;
                else
                    data_pH(i,j)=-1;
                end
            end
        end
        
        data_pH(data_pH>7.4)=7.4;
        data_pH(data_pH<6.0 & data_pH>0.0)=6.0;
end

