function [LinTblX,LinTblY]=LoadXYtr_Mex( ThetaX, ThetaY, HalfIline, HalfJline, ...
                            NumIlines, NumJlines, Dx, dx )

%tic % TICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTIC

Dy=Dx;
dy=dx;
Iimage=1024; Jimage=1024;

UpScale=64;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Létrehozzuk az uniform téglalaprács esetén minden elemben érvényes
% mátrixot, ami a polinomegyütthatókre felírt egyenletrendszert írja le.

Dxp3 = Dx / 3;
Dyp3 = Dy / 3;

% tic

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

% A mezõk tényleges feltöltése %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LinTblX LinTblY] = A40_LoadXYtr_06_Cpp( ThetaX, ThetaY, Mi, UpScale, dy, dx, HalfIline, HalfJline, Dx, Dy, NumIlines, NumJlines, Iimage, Jimage );

end