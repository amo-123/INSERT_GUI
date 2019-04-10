function [ShiftXY,Out]=A20_SplineToSnakes_08(PicX,PicY,LinesAllX1,LinesAllY1,varargin)

% tic %TICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTIC


if (nargin==6) % Vannak X2 �s Y2 vonalak is %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    LinesAllX2=varargin{1};
    LinesAllY2=varargin{2};
    
    % X vonalakat �sszef�s�lj�k %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Size1=size(LinesAllX1);
    Size2=size(LinesAllX2);
    NumJlines=Size1(1)+Size2(1);
    LinesAllX=zeros(NumJlines,1024,3);
    Odd=1;
    if(LinesAllX1(2,512,2)>LinesAllX2(2,512,2)) Odd=2; end
    for i=1:NumJlines
        if ( rem(i,2)==rem(Odd,2) )
            LinesAllX(i,:,:)=LinesAllX1((i+rem(Odd,2))/2,:,:);
        else
            LinesAllX(i,:,:)=LinesAllX2((i+1-rem(Odd,2))/2,:,:);
        end
    end
    
    % Y vonalakat �sszef�s�lj�k %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Size1=size(LinesAllY1);
    Size2=size(LinesAllY2);
    NumIlines=Size1(1)+Size2(1);
    LinesAllY=zeros(NumIlines,1024,3);
    Odd=1;
    if(LinesAllY1(2,512,2)>LinesAllY2(2,512,2)) Odd=2; end
    for i=1:NumIlines
        if ( rem(i,2)==rem(Odd,2) )
            LinesAllY(i,:,:)=LinesAllY1((i+rem(Odd,2))/2,:,:);
        else
            LinesAllY(i,:,:)=LinesAllY2((i+1-rem(Odd,2))/2,:,:);
        end
    end
    
else % Csak X1 �s Y1 vonalak vannak. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LinesAllX=LinesAllX1;
    LinesAllY=LinesAllY1;
    NumJlines=size(LinesAllX1,1);
    NumIlines=size(LinesAllY1,1);
end


% Megcsin�ljuk a t�mb�t, amibe a diszkr�t eltol�si t�r ker�l majd.
% Koordin�t�k: [Metsz�spont I][Metsz�spont J][x eltol�s, y eltol�s].
ShiftXY=zeros(NumIlines,NumJlines,2);
% Bele is t�lt�nk, ami megvan m�r vonal. Ha nem minden indexre futtatunk,
% akkor ezzel kieg�sz�thetj�k, fel�l�rhatjuk a vonal szettet.

%load('ShiftXY.mat');


SplineL=16; % Spline hossza. Lehet�leg 2 hatv�ny legyen.
SplineLp2=round(SplineL/2); SplineLp4=round(SplineL/4);


CutEnds=1;
% Ennyi pontot dobunk ki a vonalak v�g�r�l. A vonalkeres�si algoritmusban
% szerepl� meg�ll�si felt�telt�l f�gg, hogy itt mi a helyes �rt�k!

% % �t k�ne n�zni, hogy ezt hogyan is haszn�lom....
% CutJbegin=zeros(100,100);
% CutJend=zeros(100,100);
% CutIbegin=zeros(100,100);
% CutIend=zeros(100,100);
%         
%         
% for JL=1:NumJlines
%     idx=(LinesAllX(JL,:,2)>0);
%     y=[1:1024];y=y(idx);
%     if(CutJbegin(JL)~=0)
%         for i=y(1):y(CutJbegin(JL))
%             LinesAllX(JL,i,2)=0;
%         end
%     end
%     if(CutJend(JL)~=0)
%         for i=y(end-CutJend(JL)+1):y(end)
%             LinesAllX(JL,i,2)=0;
%         end
%     end
% end
% for IL=1:NumIlines
%     idx=(LinesAllY(IL,:,2)>0);
%     x=[1:1024];x=x(idx);
%     if(CutIbegin(IL)~=0)
%         for j=x(1):x(CutIbegin(IL))
%             LinesAllY(IL,j,2)=0;
%         end
%     end
%     if(CutIend(IL)~=0)
%         for j=x(end-CutIend(IL)+1):x(end)
%             LinesAllY(IL,j,2)=0;
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Megcsin�ljuk a h�tt�rk�pet, beolvassuk a felhaszn�lt felv�teleket.


% if ishandle(21) delete(21); end
% FigH=figure(21);
FigH=figure();
set(FigH,'NumberTitle','off','Name','Phantom Lines'' Crosses','MenuBar','none');
ZoomMenuH = uimenu(FigH,'Label','Zoom On/Off');
uimenu(ZoomMenuH,'Label','Zoom On','Callback','zoom on');
uimenu(ZoomMenuH,'Label','Zoom Off','Callback','zoom off');
colormap(pink);
imagesc(abs(PicX)+abs(PicY));
hold on


% Metsz�spontok keres�se %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CrossH=zeros(NumIlines,NumJlines);

for IL=1:NumIlines
    for JL=1:NumJlines
        idx=(LinesAllY(IL,:,2)>0);
        Yx=LinesAllY(IL,(idx),2);
        x=[1:1024];x=x(idx); % Csak ott vessz�k Y(x)-et, ahol �rtelmezve van.
        
        idx=(LinesAllX(JL,:,2)>0);
        Xy=LinesAllX(JL,(idx),2);
        y=[1:1024];y=y(idx);
        
        
        % Megcsin�ljuk a Spline breakpointokat. �gyel�nk, hogy a legv�g�re
        % ker�lj�n break, �s hogy az eleje �s v�ge fele/negyede hossz� legyen.
        bx=[x(1+CutEnds),x(1+SplineLp2),x(1+SplineLp2+SplineLp4)...
            :SplineL:x(end-CutEnds),0,0];
        by=[y(1+CutEnds),y(1+SplineLp2),y(1+SplineLp2+SplineLp4)...
            :SplineL:y(end-CutEnds),0,0];
        diffx=x(end-CutEnds)-bx(end-2);
        diffy=y(end-CutEnds)-by(end-2);
        for k=2:10
            bx(end-k)=bx(end-k)+diffx;
            by(end-k)=by(end-k)+diffy;
        end
        bx(end-2)=bx(end-3)+SplineLp2;
        bx(end-1)=bx(end-2)+SplineLp4;
        bx(end)=x(end-CutEnds);
        by(end-2)=by(end-3)+SplineLp2;
        by(end-1)=by(end-2)+SplineLp4;
        by(end)=y(end-CutEnds);
        
        
        plot(x,Yx,'r-')
        hold on
        plot(Xy,y,'r-')
        
        % Illeszt�nk, �br�zolunk.
        spYx=fnxtr(spline(bx,Yx/...
            spline(bx,eye(length(bx)),x)),1);
        spXy=fnxtr(spline(by,Xy/...
            spline(by,eye(length(by)),y)),1);
        plot([1:1024],fnval(spYx,[1:1024]),'g-');
        plot(fnval(spXy,[1:1024]),[1:1024],'g-');
        
        
        % Metsz�spontokat keres�nk. LorR: 'Left or Right'.
        step=512; x0=0; LorR=1; %
        x0=x0+LorR*step;
        
        while (step>0.01) % A felt�telben szerepel az �rz�kenys�g.
            y0=fnval(spYx,x0);
            xy0=fnval(spXy,y0);
            if (x0>xy0)
                if (LorR==1)  step=step/2; end
                LorR=(-1); end
            if (x0<xy0)
                if (LorR==-1) step=step/2; end
                LorR=1; end
            
            while ( ((x0+LorR*step)>1020) || ... %x(end)
                    ((x0+LorR*step<5)) )         %x(1)
                step=step/2;
            end
            x0=x0+LorR*step;
        end
        
%         CrossH(IL,JL)=plot(x0,y0,'go','MarkerSize',6); % Megjelen�tj�k a metsz�spontot k�rrel.
%         plot(xy0,y0,'g+'); % Megjelen�tj�k a metsz�spont vet�tettj�t kereszttel. Ha gond van, nem esik egybe a k�rrel.
        ShiftXY(IL,JL,1)=x0;
        ShiftXY(IL,JL,2)=y0;
        
        % Hozz�rendelj�k az x0, y0 pontokat reprezent�l� k�r�kh�z az eg�r 
        % gombnyom�sra akt�vv� tev� f�ggv�nyt.
%        set(CrossH(IL,JL),'ButtonDownFcn',{@RePosCross,IL,JL});
        
    end
end

% Defini�ljuk az x0, y0 pontokat reprezent�l� k�r�ket eg�r gombnyom�sra 
% akt�vv� tev� f�ggv�nyt.
% function RePosCross(source,eventdata,i,j)
%     set(CrossH(i,j),'Color',[1 1 0]);
%     [x,y]=getpts(gcf);
%     if (1-isempty(x))
%         set(CrossH(i,j),'XData',x(end),'YData',y(end));
%         ShiftXY(i,j,1)=x(end);
%         ShiftXY(i,j,2)=y(end);
%         save('ShiftXY.mat','ShiftXY');
%     end
%     set(CrossH(i,j),'Color',[0 1 0]);
% end

hold off

Out=toc;

end % End of function