function [ThetaX,ThetaY,Out]=SplineToCrosses(ShiftXY,HalfIline,HalfJline,...
                             NumIlines,NumJlines,Dx,dx)

% tic  %TICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTIC

Dy=Dx;
dy=dx;
Iimage=1024;
Jimage=1024;


Dxp3=Dx/3; Dyp3=Dy/3; % �lharmadol� pontok t�vols�ga.

Out=0;

DimJL=NumJlines+2; % A sz�leken extrapol�lt diszkr�t mez� m�retei.
DimIL=NumIlines+2;
DimJL3=DimJL*3-2; % A harmadol� pontok miatt 3X-os oszt�s is kell.
DimIL3=DimIL*3-2;

if (HalfIline==0)
    HalfIline=DimIL/2+0.5; % K�z�pvonalak sorsz�ma.
else
    HalfIline=HalfIline+1;
end

if (HalfJline==0)
    HalfJline=DimJL/2+0.5;
else
    HalfJline=HalfJline+1;
end

HalfIline3=HalfIline*3-2; % K�z�pvonal 3X-os oszt�sn�l.
HalfJline3=HalfJline*3-2;


% Eddig a ShiftXY pixelben adta meg a metsz�spontok koordin�t�it, innent�l
% fogva mm-ben adjuk meg a metsz�spontokban az eltol�si vektor
% komponenseket az OffsetXY t�mbben.

OffsetXY=zeros(DimIL,DimJL,2);
GridX=zeros(DimJL,1); GridY=zeros(DimIL,1);
GridX3=zeros(DimJL3+2,1); GridY3=zeros(DimIL3+2,1);
% L�trehozzuk a diszkr�t eltol�si mez�t, illetve az X �s Y r�csokat.
% GridX3 2-vel hosszabb, hogy a ciklusokn�l ne csorduljunk t�l.

for JL=1:DimJL    GridX(JL)=Dx*(JL-HalfJline)+dx*(Jimage/2-1);       end
for IL=1:DimIL    GridY(IL)=Dy*(IL-HalfIline)+dy*(Iimage/2-1);       end
for JL3=1:DimJL3  GridX3(JL3)=Dxp3*(JL3-HalfJline3)+dx*(Jimage/2-1); end
for IL3=1:DimIL3  GridY3(IL3)=Dyp3*(IL3-HalfIline3)+dy*(Iimage/2-1); end
% Ezekn�l le kell ellen�rizni, hogy a t�bbi modulhoz k�pest nem cs�sztunk-e
% meg! Jline2x �s Iline2y f�ggv�nyek a LoadXYtr modulban!!!!!!!!!!!!!!!!!!!

for IL=1:DimIL-2
    for JL=1:DimJL-2
        OffsetXY(IL+1,JL+1,1)=ShiftXY(IL,JL,1)*dx-GridX(JL+1);
        OffsetXY(IL+1,JL+1,2)=ShiftXY(IL,JL,2)*dx-GridY(IL+1);
    end
end
% Az OffsetXY indexel�s�t egyel el is kellett tolni fentebb a ShiftXY-hoz
% k�pest a sz�leken val� extrapol�ci� miatt.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extrapol�lunk a sz�leken. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extrapol�ci� "k�z�pt�jt":

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

% Extrapol�ci� a sarkok mellett:

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

% Extrapol�ci� a sarkokban:

% �tl� ment�n: Ez nagyon biztons�gosnak t�nik, de nem sz�p a sarok.
% OffsetXY(1,1,:)=2*OffsetXY(2,2,:)-OffsetXY(3,3,:);
% OffsetXY(1,DimJL,:)=2*OffsetXY(2,DimJL-1,:)-OffsetXY(3,DimJL-2,:);
% OffsetXY(DimIL,1,:)=2*OffsetXY(DimIL-1,2,:)-OffsetXY(DimIL-2,3,:);
% OffsetXY(DimIL,DimJL,:)=2*OffsetXY(DimIL-1,DimJL-1,:)-OffsetXY(DimIL-2,DimJL-2,:);

% Sz�ls� X �s Y vonalak ment�n vett extrapol�ltak �tlaga.
OffsetXY(1,1,:)=OffsetXY(2,1,:)+OffsetXY(1,2,:)...
   -0.5*OffsetXY(3,1,:)-0.5*OffsetXY(1,3,:);
OffsetXY(1,DimJL,:)=OffsetXY(2,DimJL,:)+OffsetXY(1,DimJL-1,:)...
   -0.5*OffsetXY(3,DimJL,:)-0.5*OffsetXY(1,DimJL-2,:);
OffsetXY(DimIL,1,:)=OffsetXY(DimIL-1,1,:)+OffsetXY(DimIL,2,:)...
   -0.5*OffsetXY(DimIL-2,1,:)-0.5*OffsetXY(DimIL,3,:);
OffsetXY(DimIL,DimJL,:)=OffsetXY(DimIL-1,DimJL,:)+OffsetXY(DimIL,DimJL-1,:)...
   -0.5*OffsetXY(DimIL-2,DimJL,:)-0.5*OffsetXY(DimIL,DimJL-2,:);


           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K�p l�trehoz�sa �s az OffsetXY t�mb elemeinek �br�zol�sa:

% figure(31)
% % A=zeros(round(max(NumIlines,NumJlines)*Dx*dx)+50);
% A=zeros(round(GridX(DimJL)+150),round(GridY(DimIL)+150));
% 
% imagesc(A);
% hold on
% for IL=1:DimIL
%     for JL=1:DimJL
%         plot(GridX(JL)+OffsetXY(IL,JL,1),GridY(IL)+OffsetXY(IL,JL,2),'bo');
%         if ( (IL==1) || (IL==DimIL) || (JL==1) || (JL==DimJL) )
%             plot(GridX(JL)+OffsetXY(IL,JL,1),GridY(IL)+OffsetXY(IL,JL,2),'ro');
%         end
%     end
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spline-ok l�trehoz�sa, �br�zol�sa �s Theta t�mb�k felt�lt�se %%%%%%%%%%%%

ThetaX=zeros(DimIL3+2,DimJL+2); % 2-vel hosszabb, hogy a ciklusokn�l ne csorduljunk t�l
ThetaY=zeros(DimIL3+2,DimJL+2);


% Spline-ok illeszt�se �s �br�zol�sa x tengely ment�n:

for IL=1:DimIL
    IL3=IL*3-2;
    
    spYx=spline(GridX(:),OffsetXY(IL,:,2));

    spXx=spline(GridX(:),OffsetXY(IL,:,1));

    x=[round(GridX(1)):round(GridX(DimJL))];    % �rtelmez�si tartom�ny
%     xShift=x+fnval(spXx,x);                  % Tengely menti torzul�s
%     plot(xShift,GridY(IL)+fnval(spYx,x),'r-');
    
    for JL3=1:DimJL3
        ThetaY(IL3,JL3)=fnval(spYx,GridX3(JL3));
        ThetaX(IL3,JL3)=fnval(spXx,GridX3(JL3));
    end
end

% Spline-ok illeszt�se �s �br�zol�sa y tengely ment�n:

for JL3=1:DimJL3
    spYy=spline(GridY3(1:3:DimIL3),ThetaY(1:3:DimIL3,JL3));

    spXy=spline(GridY3(1:3:DimIL3),ThetaX(1:3:DimIL3,JL3));
    
    y=[round(GridY3(1)):round(GridY3(DimIL3))];    % �rtelmez�si tartom�ny
%     yShift=y+fnval(spYy,y);                  % Tengely menti torzul�s
%     if (rem(JL3,3)==1)
%         plot(GridX3(JL3)+fnval(spXy,y),yShift,'r-');
%     end
    
    for IL3=1:DimIL3
        ThetaY(IL3,JL3)=fnval(spYy,GridY3(IL3));
        ThetaX(IL3,JL3)=fnval(spXy,GridY3(IL3));
    end
end

dThetaX=0.0;
dThetaY=0.0;
for IL3=1:DimIL3
    spYx=spline(GridX3(1:3:DimJL3),ThetaY(IL3,1:3:DimJL3));

    spXx=spline(GridX3(1:3:DimJL3),ThetaX(IL3,1:3:DimJL3));

    x=[round(GridX3(1)):round(GridX3(DimJL3))];    % �rtelmez�si tartom�ny
    
    for JL3=1:DimJL3
        ThY=fnval(spYx,GridX3(JL3));
        if ( ThetaY(IL3,JL3)~=0.0 )
            if ( abs(ThY-ThetaY(IL3,JL3)>dThetaY) )
                dThetaY=abs(ThY-ThetaY(IL3,JL3));
            end
            ThetaY(IL3,JL3)=ThY;
        end
        ThX=fnval(spXx,GridX3(JL3));
        if ( ThetaX(IL3,JL3)~=0.0 )
            if ( abs(ThX-ThetaX(IL3,JL3)>dThetaX) )
                dThetaX=abs(ThX-ThetaX(IL3,JL3));
            end
            ThetaX(IL3,JL3)=ThX;
        end
    end
end


Out=1;

%toc  %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC

end