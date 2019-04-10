function [LRF,LRFoutS]=SmoothLRF_07_par(LRFoutCell,ShiftXY,NumXlines,NumYlines,NumPMT,NumRawFiles,Filter,Dx,dx)


if (NumRawFiles==1)
    
    LRFout=LRFoutCell;
    ShiftXY=Extrapolate(ShiftXY,NumYlines+2,NumXlines+2);
    
else
    
    if (NumRawFiles==4)
        LRFoutX=LRFoutCell{1}+LRFoutCell{3};
        LRFoutY=LRFoutCell{2}+LRFoutCell{4};
    else
        LRFoutX=LRFoutCell{1};
        LRFoutY=LRFoutCell{2};
    end
    LRFout=LRFoutX+LRFoutY;
    % LRFout(:,[1:2,(NumYlines-1):NumYlines],:)=LRFoutY(:,[1:2,(NumYlines-1):NumYlines],:);
    
    for i=1:NumXlines
        for j=1:NumYlines
            for k=1:(NumPMT+1)
                if (isnan(LRFout(i,j,k))==1)
                    LRFout(i,j,k)=0;
                end
            end
        end
    end
    
    for p=1:NumPMT
        LRFout(:,:,p)=LRFout(:,:,p)./LRFout(:,:,NumPMT+1);
    end
    
    
end


LRFfineC=cell(NumPMT,1);
LRFrawC=cell(NumPMT,1);
% LRFfineSum=zeros(1024,1024);

poolobj = parpool;

parfor i=1:NumPMT
    [LRFfineC{i},LRFrawC{i}]=SmoothLRF_02(i,LRFout,ShiftXY,NumXlines,NumYlines,Filter,Dx,dx);
end
delete(poolobj);

LRF=zeros(NumPMT+1,1024,1024);
LRFoutS=zeros(NumXlines,NumYlines,NumPMT+1);

for p=1:NumPMT
    LRFfine=LRFfineC{p};
    LRF(p,:,:)=LRFfine(:,:)';
    LRFoutS(:,:,p)=LRFrawC{p};
end



for i=1:1024
    for j=1:1024
        LRF(NumPMT+1,i,j)=sum(LRF(:,i,j));
    end
end

for p=1:NumPMT
    for i=1:1024
        for j=1:1024
            LRF(p,i,j)=LRF(p,i,j)/LRF(NumPMT+1,i,j);
        end
    end
end

LRF=sqrt(LRF);

end

function OffsetXY=Extrapolate(ShiftXY,DimIL,DimJL)

OffsetXY=zeros(DimIL,DimJL,2);

for IL=1:DimIL-2
    for JL=1:DimJL-2
        OffsetXY(IL+1,JL+1,1)=ShiftXY(IL,JL,1);
        OffsetXY(IL+1,JL+1,2)=ShiftXY(IL,JL,2);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extrapolating along the sides. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for IL=3:DimIL-2
    OffsetXY(IL,1,:)=1.5*OffsetXY(IL,2,:)...
        +0.25*OffsetXY(IL-1,2,:)+0.25*OffsetXY(IL+1,2,:)...
        -0.5*OffsetXY(IL,3,:)-0.25*OffsetXY(IL-1,3,:)-0.25*OffsetXY(IL+1,3,:);
    OffsetXY(IL,DimJL,:)=1.5*OffsetXY(IL,DimJL-1,:)...
        +0.25*OffsetXY(IL-1,DimJL-1,:)+0.25*OffsetXY(IL+1,DimJL-1,:)...
        -0.5*OffsetXY(IL,DimJL-2,:)-0.25*OffsetXY(IL-1,DimJL-2,:)-0.25*OffsetXY(IL+1,DimJL-2,:);
end
for JL=3:DimJL-2
    OffsetXY(1,JL,:)=1.5*OffsetXY(2,JL,:)...
        +0.25*OffsetXY(2,JL-1,:)+0.25*OffsetXY(2,JL+1,:)...
        -0.5*OffsetXY(3,JL,:)-0.25*OffsetXY(3,JL-1,:)-0.25*OffsetXY(3,JL+1,:);
    OffsetXY(DimIL,JL,:)=1.5*OffsetXY(DimIL-1,JL,:)...
        +0.25*OffsetXY(DimIL-1,JL-1,:)+0.25*OffsetXY(DimIL-1,JL+1,:)...
        -0.5*OffsetXY(DimIL-2,JL,:)-0.25*OffsetXY(DimIL-2,JL-1,:)-0.25*OffsetXY(DimIL-2,JL+1,:);
end

% Extrapolating at the corners:

OffsetXY(2,1,:)=1.5*OffsetXY(2,2,:)+0.5*OffsetXY(3,2,:)...
               -0.5*OffsetXY(2,3,:)-0.5*OffsetXY(3,3,:);
OffsetXY(1,2,:)=1.5*OffsetXY(2,2,:)+0.5*OffsetXY(2,3,:)...
               -0.5*OffsetXY(3,2,:)-0.5*OffsetXY(3,3,:);
OffsetXY(DimIL-1,DimJL,:)=1.5*OffsetXY(DimIL-1,DimJL-1,:)...
    +0.5*OffsetXY(DimIL-2,DimJL-1,:)...
    -0.5*OffsetXY(DimIL-1,DimJL-2,:)-0.5*OffsetXY(DimIL-2,DimJL-2,:);
OffsetXY(DimIL,DimJL-1,:)=1.5*OffsetXY(DimIL-1,DimJL-1,:)...
    +0.5*OffsetXY(DimIL-1,DimJL-2,:)...
    -0.5*OffsetXY(DimIL-2,DimJL-1,:)-0.5*OffsetXY(DimIL-2,DimJL-2,:);

OffsetXY(2,DimJL,:)=1.5*OffsetXY(2,DimJL-1,:)+0.5*OffsetXY(3,DimJL-1,:)...
               -0.5*OffsetXY(2,DimJL-2,:)-0.5*OffsetXY(3,DimJL-2,:);
OffsetXY(DimIL,2,:)=1.5*OffsetXY(DimIL-1,2,:)+0.5*OffsetXY(DimIL-1,3,:)...
               -0.5*OffsetXY(DimIL-2,2,:)-0.5*OffsetXY(DimIL-2,3,:);
OffsetXY(DimIL-1,1,:)=1.5*OffsetXY(DimIL-1,2,:)+0.5*OffsetXY(DimIL-2,2,:)...
               -0.5*OffsetXY(DimIL-1,3,:)-0.5*OffsetXY(DimIL-2,3,:);
OffsetXY(1,DimJL-1,:)=1.5*OffsetXY(2,DimJL-1,:)+0.5*OffsetXY(2,DimJL-2,:)...
               -0.5*OffsetXY(3,DimJL-1,:)-0.5*OffsetXY(3,DimJL-2,:);


OffsetXY(1,1,:)=OffsetXY(2,1,:)+OffsetXY(1,2,:)...
   -0.5*OffsetXY(3,1,:)-0.5*OffsetXY(1,3,:);
OffsetXY(1,DimJL,:)=OffsetXY(2,DimJL,:)+OffsetXY(1,DimJL-1,:)...
   -0.5*OffsetXY(3,DimJL,:)-0.5*OffsetXY(1,DimJL-2,:);
OffsetXY(DimIL,1,:)=OffsetXY(DimIL-1,1,:)+OffsetXY(DimIL,2,:)...
   -0.5*OffsetXY(DimIL-2,1,:)-0.5*OffsetXY(DimIL,3,:);
OffsetXY(DimIL,DimJL,:)=OffsetXY(DimIL-1,DimJL,:)+OffsetXY(DimIL,DimJL-1,:)...
   -0.5*OffsetXY(DimIL-2,DimJL,:)-0.5*OffsetXY(DimIL,DimJL-2,:);

end
