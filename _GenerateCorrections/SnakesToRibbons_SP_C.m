function [LinesAllMat,NumLines]=SnakesToRibbons_SP_C(Pic,Orient,...
              MidLine,Extend,LeftSide,RightSide,FindOnly,Tight)

% Built-in parameters. To be genralized later.
NumIlines=100; NumJlines=100;
Iimage=1024; Jimage=1024;
Sigma=5;


LinesAllMat=0;

if (Orient=='X')
    NumLines=NumJlines;
else
    NumLines=NumIlines;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating background image.

Pic=abs(Pic);
if ( (strcmp(Orient,'Y')==1) ) Pic=Pic'; end

colormap(pink);
imagesc(Pic);
hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start seed search. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function PeaksAtMidLine=LineFinder()
    
    if ( (strcmp(Orient,'X')==1) ) 
        MinimumPeak = sum(sum(Pic(212:811,362:661)))/600/300 ;
    else
        MinimumPeak = sum(sum(Pic(362:661,212:811)))/600/300 ;
    end
    
    MidVector=sum(Pic((MidLine-5):(MidLine+5),:)) / 9;
    
    %clinical trying
%     filter=[0.2:0.2:0.8,1,0.8:-0.2:0.2]; %/sum([0.2:0.2:0.8,1,0.8:-0.2:0.2]);
    filter=[0.4:0.3:1,0.7:-0.3:0.4]; %/sum([0.4:0.3:1,0.7:-0.3:0.4]);

    %preclinical opt:
%     filter=[0.1:0.1:0.9,1,0.9:-0.1:0.1]; % Original filter
    
    %old comment:
    % filter=[0,0,0,0,0,0.2:0.2:0.8,1,0.8:-0.2:0.2,0,0,0,0,0];
    % filter=[0,0,0,0,0,0,0,0.4:0.3:0.7,1,0.7:-0.3:0.4,0,0,0,0,0,0,0];
    mask=zeros(1,5);
    FilteredValue=zeros(1,1024);
    PeaksAtMidLine=zeros(1,NumLines);

    Peak=1;
    if ( strcmp(Extend,'Left') || strcmp(Extend,'Both') )
        PeaksAtMidLine(1)=LeftSide;
        Peak=2;
    end
    
    left=20;
    if((strcmp(Extend,'Left')==1)||(strcmp(Extend,'Both')==1))
        left=LeftSide+5;
        FilteredValue(left-2)=1e5; FilteredValue(left-1)=0.9e5;
    end
    right=1005;
    if((strcmp(Extend,'Right')==1)||(strcmp(Extend,'Both')==1))
        right=RightSide-5;
    end
    
    jump=0;
    for j=left:right
%         mask(:)=MidVector( (j-9):(j+9) )'; preclinical
        mask(:)=MidVector( (j-2):(j+2) )';
        FilteredValue(j)=filter*mask';
        if jump==0
            if( (FilteredValue(j)<FilteredValue(j-1))...
                    && (FilteredValue(j-1)>=FilteredValue(j-2))...
                    && (MidVector(j) > MinimumPeak ) )
                PeaksAtMidLine(Peak)=j-1;
                Peak=Peak+1;
                jump=2;
            end
        else
            jump=jump-1;
        end
    end
    
    NumLines=Peak+strcmp(Extend,'Right')+strcmp(Extend,'Both')-1;
    if ( strcmp(Extend,'Right') || strcmp(Extend,'Both') )
        PeaksAtMidLine(NumLines)=RightSide;
    end
    
    PeaksAtMidLine=PeaksAtMidLine(1:NumLines);
    
    if (NumLines>41); 
        PeaksAtMidLine=PeaksAtMidLine(2:end); 
        NumLines=NumLines-1; 
    end  
    
    plot(PeaksAtMidLine,MidLine,'yo')
    
end


% Calling start seed search. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PeaksAtMidLine=LineFinder();

% Return, if it is in FindOnly mode. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (strcmp(FindOnly,'TRUE'))
    return
end


% On which side to be used Tight Search?

if ( strcmp(Tight,'None'))  TightL=0; TightR=0; end
if ( strcmp(Tight,'Left'))  TightL=1; TightR=0; end
if ( strcmp(Tight,'Right')) TightL=0; TightR=1; end
if ( strcmp(Tight,'Both'))  TightL=1; TightR=1; end


% End of initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LinesAllMat=zeros(NumLines,1024,3);
LinesAllCell=num2cell(LinesAllMat);
LinesAllMatM=zeros(NumLines,1024,3);
LinesAllCellM=num2cell(LinesAllMat);
LinesAllMatP=zeros(NumLines,1024,3);
LinesAllCellP=num2cell(LinesAllMat);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main line search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J1=2; Jend=NumLines-1; %J1=2; Jend=NumLines-1;
% % % % tempPic=Pic;
% % % % pf=figure;
% % % % imagesc(Pic);
% % % % darab=round(input('darab:'));
% % % % for idb=1:1:darab;
% % % % wwwpar=round(ginput(2));
% % % % for iwww=wwwpar(1,2):wwwpar(2,2)
% % % %     Pic(iwww,wwwpar(1,1):wwwpar(2,1)) = Pic(wwwpar(1,2),wwwpar(1,1):wwwpar(2,1)) ;
% % % % end
% % % % end
% % % % imagesc(Pic);
% % % % pause
% % % % close(pf)

% Pic and PeaksAtMidLine is written out for the exe.
dlmwrite('_GenerateCorrections/image.raw',Pic,' ');
dlmwrite('_GenerateCorrections/peaks.txt',PeaksAtMidLine,' ');
% % % % Pic=tempPic;
% Stop criteria.

    if ( (strcmp(Orient,'X')==1) ) 
        StopFinder = sum(sum(Pic(212:811,362:661)))/600/300 /4;
    else
        StopFinder = sum(sum(Pic(362:661,212:811)))/600/300 /4;
    end

cmdstr = ['START /WAIT ' '_GenerateCorrections/LinTabGen3GaussFit.exe ' num2str(MidLine) ' '  num2str(Sigma) ...
    ' ' num2str(J1) ' '  num2str(Jend) ' ' num2str(NumLines) ' ' '_GenerateCorrections/peaks.txt' ...
    ' ' num2str(1024) ' ' num2str(1024) ' ' '_GenerateCorrections/image.raw' ' ' num2str(TightL) ...
    ' ' num2str(TightR) ' '  num2str(StopFinder) ' ' '1e-5 ' '40 ' '1000' ...
    ' '  num2str(1) ' ' '_GenerateCorrections/' ]; % ' '  num2str(3) 

disp(cmdstr);

s = system(cmdstr);

% Middle lines: Lines; Minus (left) side lines: LinesM; Positive (right)
% side lines: LinesP;
Lines  = dlmread( '_GenerateCorrections/lines.txt' );
LinesM = dlmread( '_GenerateCorrections/linesM.txt');
LinesP = dlmread( '_GenerateCorrections/linesP.txt');

LinesAllMat=zeros(NumLines,1024,3); LinesAllMatM=LinesAllMat; LinesAllMatP=LinesAllMat;

for JColumn = J1:Jend
    LinesAllMat (JColumn,:,:) = Lines (:,(1+(JColumn-1)*3):(3+(JColumn-1)*3));
    LinesAllMatM(JColumn,:,:) = LinesM(:,(1+(JColumn-1)*3):(3+(JColumn-1)*3));
    LinesAllMatP(JColumn,:,:) = LinesP(:,(1+(JColumn-1)*3):(3+(JColumn-1)*3));
end

LinesAllMat((J1-1),:,:)=LinesAllMatM(J1,:,:);
LinesAllMat((Jend+1),:,:)=LinesAllMatP(Jend,:,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of main line search %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


LinesAllMat2=LinesAllMat;

PutLines2Screen();


% Putting the lines into the figure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function PutLines2Screen()
        Pic2=Pic;
        for Jcolumn=1:NumLines
            for Irow=1:1024; % I1:step:Iend
                b2=LinesAllMat2(Jcolumn,Irow,2);
                if (b2>0)
                    Pic2(Irow,round(b2))=Pic2(Irow,round(b2))+100; %b2
                end
            end
        end
        imagesc(Pic2);
        axis([212 811 212 811]);
        for Jcolumn=1:NumLines
            DotH(Jcolumn)=plot(PeaksAtMidLine(Jcolumn),MidLine,'yo');
            set(DotH(Jcolumn),'ButtonDownFcn',{@RePosLine,Jcolumn});
        end
    end

end
