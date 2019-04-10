function LUT=LUT_Generator_02(ShiftXY, Rad, NumRawFiles)

dims=size(ShiftXY);

LUTX1=zeros(1024,1024,2);
LUTX2=LUTX1;
LUTY1=LUTX1;
LUTY2=LUTX1;

modulo=mod((dims-1)/2,2);

for LUTi = 1:dims(1)
    for LUTj = 1:dims(2)
        
        i=round(ShiftXY(LUTi,LUTj,2));
        j=round(ShiftXY(LUTi,LUTj,1));
        for ri=(-Rad):Rad
            for rj=(-Rad):Rad
                
                if ((ri^2+rj^2)<=Rad^2)
                    if ( (mod(LUTj,2)==modulo(2)) || (NumRawFiles==2) )
                        LUTX1(i+ri,j+rj,1)=LUTi;
                        LUTX1(i+ri,j+rj,2)=LUTj;
                    else
                        LUTX2(i+ri,j+rj,1)=LUTi;
                        LUTX2(i+ri,j+rj,2)=LUTj;
                    end
                    
                    if ( (mod(LUTi,2)==modulo(1)) || (NumRawFiles==2) )
                        LUTY1(i+ri,j+rj,1)=LUTi;
                        LUTY1(i+ri,j+rj,2)=LUTj;
                    else
                        LUTY2(i+ri,j+rj,1)=LUTi;
                        LUTY2(i+ri,j+rj,2)=LUTj;
                    end
                    
                end
            end
        end
    end
end

LUT{1}=LUTX1;
LUT{2}=LUTY1;
if (NumRawFiles==4)
    LUT{3}=LUTX2;
    LUT{4}=LUTY2;
end

% % str='LUTX1';
% % save(str,'LUTX1');
% % str='LUTX2';
% % save(str,'LUTX2');
% figure();
% Pic=LUTX1(:,:,1);
% imagesc(Pic);
% figure();
% Pic=LUTX1(:,:,2);
% imagesc(Pic);
% figure();
% Pic=LUTX2(:,:,1);
% imagesc(Pic);
% figure();
% Pic=LUTX2(:,:,2);
% imagesc(Pic);
% 
% % str='LUTY1';
% % save(str,'LUTY1');
% % str='LUTY2';
% % save(str,'LUTY2');
% figure();
% Pic=LUTY1(:,:,1);
% imagesc(Pic);
% figure();
% Pic=LUTY1(:,:,2);
% imagesc(Pic);
% figure();
% Pic=LUTY2(:,:,1);
% imagesc(Pic);
% figure();
% Pic=LUTY2(:,:,2);
% imagesc(Pic);


end