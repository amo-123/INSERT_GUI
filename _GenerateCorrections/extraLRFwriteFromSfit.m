function [LRF,PMTxy]=extraLRFwriteFromSfit(mmPerPix,NumPM,path)

% Loading a mat file which contains LRF parameters comming from Michele's
% method.
% disp('LRFwriteFromSfit :: Load LRF data...')


load(path)

% disp('LRFwriteFromSfit :: Creating starting LRF set...')



LRF=zeros(NumPM+1,1024,1024);



PMTxy=zeros(NumPM,2);



for p=1:NumPM
    
    a = LRFs{p}.a;
    b = LRFs{p}.b;
    c = LRFs{p}.c;
    offset = LRFs{p}.offset;
    x0 = LRFs{p}.x0;
    y0 = LRFs{p}.y0;

    PMTxy(p,1) = (x0-25) / mmPerPix + 512;
    PMTxy(p,2) = (y0-25) / mmPerPix + 512;
    
    for i=1:1024
        
       xmx0 = (i - 512) * mmPerPix - (x0-25);
        
        for j=1:1024
            
            ymy0 = (j - 512) * mmPerPix - (y0-25);
                        
            LRF(p,i,j) = a*exp(-(b*(xmx0)^2+c*(ymy0)^2)) + offset;
            
        end
    end
    
end



% Got the sqrt of LRF before normalization.
LRF=sqrt(LRF);

% Normalization
for i=1:1024
    for j=1:1024
        LRF(NumPM+1,i,j)=sum(LRF(:,i,j));
    end
end

for p=1:NumPM
    for i=1:1024
        for j=1:1024
            LRF(p,i,j)=LRF(p,i,j)/LRF(NumPM+1,i,j);
        end
    end
end


% Get the sqrt again after normalization.
LRF=sqrt(LRF);

end