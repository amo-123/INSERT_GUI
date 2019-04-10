function [ShiftXY,Out]=A20_SplineToSnakes_08(PicX,PicY,LinesAllX1,LinesAllY1,varargin)

% tic %TICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTIC


if (nargin==6) % Vannak X2 és Y2 vonalak is %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    LinesAllX2=varargin{1};
    LinesAllY2=varargin{2};
    
    % X vonalakat összefésüljük %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    % Y vonalakat összefésüljük %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
else % Csak X1 és Y1 vonalak vannak. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LinesAllX=LinesAllX1;
    LinesAllY=LinesAllY1;
    NumJlines=size(LinesAllX1,1);
    NumIlines=size(LinesAllY1,1);
end


% Megcsináljuk a tömböt, amibe a diszkrét eltolási tér kerül majd.
% Koordináták: [Metszéspont I][Metszéspont J][x eltolás, y eltolás].
ShiftXY=zeros(NumIlines,NumJlines,2);
% Bele is töltünk, ami megvan már vonal. Ha nem minden indexre futtatunk,
% akkor ezzel kiegészíthetjük, felülírhatjuk a vonal szettet.

%load('ShiftXY.mat');


SplineL=16; % Spline hossza. Lehetõleg 2 hatvány legyen.
SplineLp2=round(SplineL/2); SplineLp4=round(SplineL/4);


CutEnds=1;
% Ennyi pontot dobunk ki a vonalak végérõl. A vonalkeresési algoritmusban
% szereplõ megállási feltételtõl függ, hogy itt mi a helyes érték!

% % Át kéne nézni, hogy ezt hogyan is használom....
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
% Megcsináljuk a háttérképet, beolvassuk a felhasznált felvételeket.


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


% Metszéspontok keresése %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CrossH=zeros(NumIlines,NumJlines);

for IL=1:NumIlines
    for JL=1:NumJlines
        idx=(LinesAllY(IL,:,2)>0);
        Yx=LinesAllY(IL,(idx),2);
        x=[1:1024];x=x(idx); % Csak ott vesszük Y(x)-et, ahol értelmezve van.
        
        idx=(LinesAllX(JL,:,2)>0);
        Xy=LinesAllX(JL,(idx),2);
        y=[1:1024];y=y(idx);
        
        
        % Megcsináljuk a Spline breakpointokat. Ügyelünk, hogy a legvégére
        % kerüljön break, és hogy az eleje és vége fele/negyede hosszú legyen.
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
        
        % Illesztünk, ábrázolunk.
        spYx=fnxtr(spline(bx,Yx/...
            spline(bx,eye(length(bx)),x)),1);
        spXy=fnxtr(spline(by,Xy/...
            spline(by,eye(length(by)),y)),1);
        plot([1:1024],fnval(spYx,[1:1024]),'g-');
        plot(fnval(spXy,[1:1024]),[1:1024],'g-');
        
        
        % Metszéspontokat keresünk. LorR: 'Left or Right'.
        step=512; x0=0; LorR=1; %
        x0=x0+LorR*step;
        
        while (step>0.01) % A feltételben szerepel az érzékenység.
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
        
%         CrossH(IL,JL)=plot(x0,y0,'go','MarkerSize',6); % Megjelenítjük a metszéspontot körrel.
%         plot(xy0,y0,'g+'); % Megjelenítjük a metszéspont vetítettjét kereszttel. Ha gond van, nem esik egybe a körrel.
        ShiftXY(IL,JL,1)=x0;
        ShiftXY(IL,JL,2)=y0;
        
        % Hozzárendeljük az x0, y0 pontokat reprezentáló körökhöz az egér 
        % gombnyomásra aktívvá tevõ függvényt.
%        set(CrossH(IL,JL),'ButtonDownFcn',{@RePosCross,IL,JL});
        
    end
end

% Definiáljuk az x0, y0 pontokat reprezentáló köröket egér gombnyomásra 
% aktívvá tevõ függvényt.
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