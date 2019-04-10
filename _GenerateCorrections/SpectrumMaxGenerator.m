function SpectrumMax=SpectrumMaxGenerator(SpectrumCell,NumXlines,NumYlines,NumRawFiles)

if (NumRawFiles==4)
    SpectrumX=SpectrumCell{1}+SpectrumCell{3};
    SpectrumY=SpectrumCell{2}+SpectrumCell{4};
else
    SpectrumX=SpectrumCell{1};
    SpectrumY=SpectrumCell{2};
end

% Spectrum=zeros(NumXlines,NumYlines,500);
% Spectrum(5:(NumXlines-4),:,:)=SpectrumX(5:(NumXlines-4),:,:);
% Spectrum([1:4,(NumXlines-3):NumXlines],:,:)=SpectrumY([1:4,(NumXlines-3):NumXlines],:,:);

Spectrum=SpectrumX+SpectrumY;

SpectrumRatio=zeros(NumXlines,NumYlines);

for i=1:NumXlines
   for j=1:NumYlines
       SpectrumRatio=sum(SpectrumX(i,j,:))/sum(SpectrumY(i,j,:));
       if SpectrumRatio>3
           Spectrum(i,j,:)=SpectrumX(i,j,:);
       end
       if SpectrumRatio<0.333
           Spectrum(i,j,:)=SpectrumY(i,j,:);
       end
   end
end


for i=1:NumXlines
   for j=1:NumYlines
       for k=1:500
           if (isnan(Spectrum(i,j,k))==1)
               Spectrum(i,j,k)=0;
           end
       end
   end
end


SpPl=1:500;
Ch=1:500;

SpectrumMax=zeros(NumXlines,NumYlines);

for i=1:NumXlines
    for j=1:NumYlines
        SpPl(:)=Spectrum(i,j,:);
        MaxVal=max(SpPl(15:300));
        Idx=(SpPl==MaxVal);
        if (sum(Idx)>0) && sum(SpPl(15:300))>0
            MaxAt=Ch(Idx);
            % MaxAt=sum(MaxAt)/sum(Idx);
            MaxAt=max(MaxAt);
            SpectrumMax(i,j)=MaxAt*100;
%             if (sum(Idx)>1)
%             disp(['i=',num2str(i),' j=',num2str(j),...
%                 ' MaxVal=',num2str(MaxVal),' MaxAt=',num2str(MaxAt)]);
%             end
        end
    end
end

% return;

% figure()
% imagesc(SpectrumMax);


%%%%%%%%
Gauss='a1*exp(-((x-b1)/c1)^2)';
f=fittype(Gauss);


for i=1:NumXlines
    for j=1:NumYlines
        X=[ round(SpectrumMax(i,j)/100*0.8):round(SpectrumMax(i,j)/100*1.2) ];
        SpPl(:)=Spectrum(i,j,:);
        s=fitoptions('Method','NonlinearLeastSquares',...
            'StartPoint',[100,SpectrumMax(i,j)/100,10]);
        try
            gfit=fit(X',SpPl(X)',f,s);
            SpectrumMax(i,j)=gfit.b1*100;
        catch exception
            disp([num2str(i),' ',num2str(j)]);
            SpectrumMax(i,j)=0;
        end
        
        
    end
end
%%%%%%%%
end

% save('SpectrumMax','SpectrumMax');

% SpectrumMax(10,10)
% 
% 
% a1=gfit.a1
% b1=gfit.b1
% c1=gfit.c1
% 
% Y=a1*exp(-((X-b1)/c1).^2);
% 
% figure(); % Spektrum ábrázolása %%%%%%%%%%%%%%
% plot([1:500],SpPl,'b+');
% hold on
% 
% plot(1:500,Y,'r-');