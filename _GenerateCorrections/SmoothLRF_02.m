function [LRFfine,LRFraw]=SmoothLRF_02(IndPMT,LRFout,ShiftXY,NumXlines,NumYlines,Filter,Dx,dx)


Dy=Dx;
dy=dx;
DimIL=NumXlines+2;
DimJL=NumYlines+2;
Iimage=1024;
Jimage=Iimage;

UpScale=1;
IimageUp=Iimage*UpScale; JimageUp=Jimage*UpScale;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dxp3=Dx/3; Dyp3=Dy/3; % Élharmadoló pontok távolsága.

DimJL3=DimJL*3-2; % A harmadoló pontok miatt 3X-os osztás is kell.
DimIL3=DimIL*3-2;

HalfIline=DimIL/2+0.5; % Középvonalak sorszáma.
HalfJline=DimJL/2+0.5;
HalfIline3=HalfIline*3-2; % Középvonal 3X-os osztásnál.
HalfJline3=HalfJline*3-2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LRFp=zeros(NumXlines+2,NumYlines+2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zajszûrés !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


% Kisegítõ függvények %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Iline=Ipix2Iline(Ipix)
        Iline = dy * ( (Ipix-0.5) - (Iimage/2-1) ) / Dy + HalfIline - 1;
    end

    function Jline=Jpix2Jline(Jpix)
        Jline = dx * ( (Jpix-0.5) - (Jimage/2-1) ) / Dx + HalfJline - 1;
    end

% End of kisegít? függvények %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:NumXlines
    for j=1:NumYlines
        Xp(i+NumXlines*(j-1))=i;
        Yp(i+NumXlines*(j-1))=j;
        Zp(i+NumXlines*(j-1))=LRFout(i,j,IndPMT);
    end
end
    

[xData, yData, zData] = prepareSurfaceData( Xp, Yp, Zp );

ft = fittype( 'loess' );
opts = fitoptions( ft );
opts.Span = Filter;
opts.Normalize = 'on';

[fitresult, gof] = fit( [xData, yData], zData, ft, opts );

if ShiftXY==0;
    for i=1:DimIL
        for j=1:DimJL
            LRFp(i,j)=fitresult(i-1,j-1);
        end
    end
else
    for i=1:DimIL
        for j=1:DimJL
            imod=Ipix2Iline(ShiftXY(j,i,1));
            jmod=Jpix2Jline(ShiftXY(j,i,2));
            LRFp(i,j)=fitresult(imod,jmod);
        end
    end
end


for IL=1:DimIL
    for JL=1:DimJL
        if (isnan(LRFp(IL,JL))) LRFp(IL,JL)=0; end
%         if (LRFp(IL,JL,NumPMT+1)<1000) LRFp(IL,JL,IndPMT)=0; end % Kinullázzuk, ahol kevés a beütés
    end
end


% Zajszûrés !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LRFraw=LRFp( 2:(DimIL-1) , 2:(DimJL-1) );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




GridX=zeros(DimJL,1); GridY=zeros(DimIL,1);
GridX3=zeros(DimJL3+2,1); GridY3=zeros(DimIL3+2,1);
% Létrehozzuk a diszkrét eltolási mezõt, illetve az X és Y rácsokat.
% GridX3 2-vel hosszabb, hogy a ciklusoknál ne csorduljunk túl.

for JL=1:DimJL    GridX(JL)=Dx*(JL-HalfJline)+dx*(Jimage/2-1);       end
for IL=1:DimIL    GridY(IL)=Dy*(IL-HalfIline)+dy*(Iimage/2-1);       end
for JL3=1:DimJL3  GridX3(JL3)=Dxp3*(JL3-HalfJline3)+dx*(Jimage/2-1); end
for IL3=1:DimIL3  GridY3(IL3)=Dyp3*(IL3-HalfIline3)+dy*(Iimage/2-1); end
           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spline-ok létrehozása, ábrázolása és a Theta tömb feltöltése %%%%%%%%%%%%

Theta=zeros(DimIL3+2,DimJL+2); % 2-vel hosszabb, hogy a ciklusoknál ne csorduljunk túl

% Spline-ok illesztése és ábrázolása x tengely mentén:

for IL=1:DimIL
    IL3=IL*3-2;
    
    Vect=1:DimJL;
    Idx=(LRFp(IL,Vect)>0);
    Vect=Vect(Idx);
    spThetaX=spline(GridX(Vect),LRFp(IL,Vect));

    for JL3=1:DimJL3
        Theta(IL3,JL3)=fnval(spThetaX,GridX3(JL3));
    end
end

% Spline-ok illesztése és ábrázolása y tengely mentén:

for JL3=1:DimJL3
    spThetaY=spline(GridY3(1:3:DimIL3),Theta(1:3:DimIL3,JL3));
    
    for IL3=1:DimIL3
        Theta(IL3,JL3)=fnval(spThetaY,GridY3(IL3));
    end
end

dThetaX=0.0;
for IL3=1:DimIL3
    spThetaX=spline(GridX3(1:3:DimJL3),Theta(IL3,1:3:DimJL3));

    
    for JL3=1:DimJL3
        ThX=fnval(spThetaX,GridX3(JL3));
        if ( Theta(IL3,JL3)~=0.0 )
            if ( abs(ThX-Theta(IL3,JL3)>dThetaX) )
                dThetaX=abs(ThX-Theta(IL3,JL3));
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A40_LoadXYtr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Létrehozzuk az uniform téglalaprács esetén minden elemben érvényes
% mátrixot, ami a polinomegyütthatókre felírt egyenletrendszert írja le.

M=zeros(16,16);
M(:,1)=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]';
M(:,2)=Dxp3*[0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3]';
M(:,3)=Dyp3*[0,0,0,0,1,1,1,1,2,2,2,2,3,3,3,3]';
M(:,4)=M(:,2).^2;
M(:,5)=M(:,2).*M(:,3);
M(:,6)=M(:,3).^2;
M(:,7)=M(:,2).^3;
M(:,8)=M(:,2).^2.*M(:,3);
M(:,9)=M(:,2).*M(:,3).^2;
M(:,10)=M(:,3).^3;
M(:,11)=M(:,2).^3.*M(:,3);
M(:,12)=M(:,2).^2.*M(:,3).^2;
M(:,13)=M(:,2).*M(:,3).^3;
M(:,14)=M(:,2).^3.*M(:,3).^2;
M(:,15)=M(:,2).^2.*M(:,3).^3;
M(:,16)=M(:,2).^3.*M(:,3).^3;

% Létrehozzuk az inverz mátrixot.
Mi=M^-1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Létrehozzuk és feltöltjük az LRFfine mezõt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LRFfine=zeros(IimageUp,JimageUp);

for IL=1:(DimIL-1)
    IL3=IL*3-2;
    Ibegin=ceil( ( Dy*(IL-HalfIline)+dy*(Iimage/2-1) )*UpScale/dy );
    Iend=ceil( ( Dy*((IL+1)-HalfIline)+dy*(Iimage/2-1) )*UpScale/dy ) -1;
    
%     Vect=1:DimJL;
%     Idx=(LRFp(IL,Vect,IndPMT)>0);
%     Vect=Vect(Idx);
%     Vect=(min(Vect)):(max(Vect)-1);
   
    for JL=1:(DimJL-1) % JL=Vect  % JL=1:(DimJL-1)
        JL3=JL*3-2;
        Jbegin=ceil( ( Dx*(JL-HalfJline)+dx*(Jimage/2-1) )*UpScale/dx );
        Jend=ceil( ( Dx*((JL+1)-HalfJline)+dx*(Jimage/2-1) )*UpScale/dx ) -1;

        % A külsõ ciklusokat a lineáris fantom vonalak által kijelölt
        % elemek szintjén szervezzük meg.
        
        % Kiszámítjuk az adott elemben az együttható vektort.
        
        T=zeros(1,16);
        T(1:4) = Theta(IL3,  JL3:(JL3+3));
        T(5:8) = Theta(IL3+1,JL3:(JL3+3));
        T(9:12)= Theta(IL3+2,JL3:(JL3+3));
        T(13:16)=Theta(IL3+3,JL3:(JL3+3));
        alfa=Mi*T';

        % Az elemen belül végimegyünk minden pixelen és feltöltjük az LRF mezõt
        
        for I=Ibegin:Iend
            for J=Jbegin:Jend
                
                % x és y értéke az elem sarkától számítva. Ez a feltétele
                % annak, hogy egy mátrixszal írjuk le az összes elemet.
                x=dx*(J-0.5)/UpScale-(Dx*(JL-HalfJline)+dx*(Jimage/2-1));
                y=dy*(I-0.5)/UpScale-(Dy*(IL-HalfIline)+dy*(Iimage/2-1));
           
                LRFfine(I,J)=...
                alfa(1)+alfa(2)*x+alfa(3)*y+alfa(4)*x^2+...
                alfa(5)*x*y+alfa(6)*y^2+alfa(7)*x^3+alfa(8)*x^2*y+...
                alfa(9)*x*y^2+alfa(10)*y^3+alfa(11)*x^3*y+...
                alfa(12)*x^2*y^2+alfa(13)*x*y^3+...
                alfa(14)*x^3*y^2+alfa(15)*x^2*y^3+alfa(16)*x^3*y^3;
            end
        end
    end
    
end

% figure();
% imagesc(LRFfine);

end