function histogram(Image,nbins,min,max,display)

for i=1:1:size(InputIm,1)
	     for j=1:1:size(InputIm,2)      
		value=InputIm(i,j); %read input image level    
		if value >= minvalue && value <= maxvalue
		    for i=1:1:nbins 
		        if value <= i*binsize+minvalue && value > (i-1)*binsize+minvalue   
			    % original histogram in pixels
		            InputIm_histogram(i)=InputIm_histogram(i)+1;
			    % normalized histogram pdf  
		            InputIm_normalized_histogram(i)=InputIm_histogram(i)/resolution;   
		        end
		    end
		end
	    end
end