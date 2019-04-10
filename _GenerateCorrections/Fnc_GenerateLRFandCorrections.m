function Fnc_GenerateLRFandCorrections(Head, NumXlines, NumYlines, Dx, dx, ...
             NumPM, pathXfile, pathYfile, pathFloodFile, WorkDir)

for iHead=Head:Head %1:10
             Head=iHead;
             pathXfile = [WorkDir '\Corr\digitalraw\digitalraw_H' sprintf('%0.2d',Head) '_X.dat'];
             pathYfile = [WorkDir '\Corr\digitalraw\digitalraw_H' sprintf('%0.2d',Head) '_Y.dat'];
             pathFloodFile = [WorkDir '\Corr\digitalraw\digitalraw_H' sprintf('%0.2d',Head) '_F.dat'];
             

tic; %TICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTIC
disp('HEREEEEEE');

NumRawFiles=2; % It means that there is <X+Y> acquisitions. Originally, this script was created for both <X+Y> and <X1+X2+Y1+Y2> acquisitions.


% try
%     load _TMP\MidLine.mat % See documentation
%     MidLine1;
%     MidLine2;
%     disp(['Midlines are set by MidLine.mat: MidLine1 = ',num2str(MidLine1), ', Midline2 = ', num2str(MidLine2),'.']);
% catch
%     MidLine1 = 512;
%     MidLine2 = 512;
%     
% end

if dx==0.2
    MidLine1 = 512;
    MidLine2 = 512;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main parameters
disp('modified in an random way')
BaseLine_vektor=[113 105 105 105 105 115 115 85 115 115 110 110 110 110 110 110 110 110 110 110]; % These are the optimized values for calibration acq at 2017-10-02
BaseLine=BaseLine_vektor(Head);

Filter1=0.2; % Filtering of the LRFs generated based on measurement data.gui
% Filter2=0.01; % Removed from original script.

LUT_mask1=10; % Size of the dots on LUT table.
% LUT_mask2=3; % Removed from original script.

% LRFtrying

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Starting LRF-s from Michele's algorithm, and PMT coordinates from the
% same MatLab file.
disp('HERE')
[LRF,PMTxy]=LRFwriteFromSfit(dx,NumPM,Head);
% LRFstart=LRF;
% if dx==0.2
%     load '_TMP\LRFmatrixs.mat' %1024*1024-es, de 0.2mm pixel, module1-hez adott olasz LRF-vel
% else
%     load '_TMP\LRFmatrixs_double.mat' %2048*1024-as mátrixhoz, 0.1mm, module1-hez adott olasz LRF-vel
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data files

GeneralLoop=0; % If it is 0, then we go through the whole file.

path{1}=pathXfile;
path{2}=pathYfile;
for PhantomIndex=1:NumRawFiles
    
    StreamFile=fopen(path{PhantomIndex},'r');
    fseek(StreamFile,0,'eof');
    Loop(PhantomIndex)=floor(ftell(StreamFile)/(72*2+4*2));
    fclose(StreamFile);
    
    % If GeneralLoop is 0, then we go through the whole file.
    if ( GeneralLoop~=0 && GeneralLoop<Loop(PhantomIndex) )
        Loop(PhantomIndex)=GeneralLoop;
    end
    
    % Rounding by 1e4. See StreamDataLength values in MEX C files.
    Loop(PhantomIndex)=Loop(PhantomIndex)-rem(Loop(PhantomIndex),1e4);
    
    % Path is rewritten for C
    path{PhantomIndex} = regexprep(path{PhantomIndex},'\','/');

end

StreamFile=fopen(pathFloodFile,'r');
fseek(StreamFile,0,'eof');
LoopUC=floor(ftell(StreamFile)/(72*2+4*2));
fclose(StreamFile);

% Rounding by 1e4.
LoopUC=LoopUC-rem(LoopUC,1e4);

% Path rewritten for C.
pathFloodFile = regexprep(pathFloodFile,'\','/');


% Initialization of main cell variables.
SpectrumCell=cell(NumRawFiles,1);
LRFoutCell=cell(NumRawFiles,1);
PicCell=cell(NumRawFiles,1);
LUT=cell(NumRawFiles,1);
for PhantomIndex=1:NumRawFiles
    LUT{PhantomIndex}=zeros(1024,1024,2);
end

SpectrumMax=zeros(NumXlines,NumYlines); % Nincs jelentõsége, nem fogjuk használni, csak kell input a MEX állománynak.
PE=17000;

% FIRST PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acquisitions are run with Statistic Method without energy windowing.
% Then, the LUT tables are created with the results.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc; %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC
disp(' ');
disp('NULL PART: Playing acquisitions with Anger (centroid) Method for looking at the pre-images (BaseLine checking).');
% Play the acquisition %%%%%%%%%%%%%%%%%%%%%
SpectrumWindow=-1; % Which means that there is no energy window and LUTs are not used.
for PhantomIndex=1:NumRawFiles
[PicCell{PhantomIndex},Count(PhantomIndex),CountEw,LRFoutCell{PhantomIndex},SpectrumCell{PhantomIndex}]=...
    MexSPEngine_08insert_ANGER(LRF, path{PhantomIndex}, Loop(PhantomIndex), PMTxy, SpectrumMax,...
    SpectrumWindow, LUT{PhantomIndex}, NumXlines, NumYlines, NumPM,PE,BaseLine);
disp(['Count: ',num2str(Count(PhantomIndex))]);
end
% Figure for checking the anger
FigH=figure();
set(FigH,'Position',[0 0 1920 1080]);
set(FigH,'NumberTitle','off','Name','Anger images');
for i=1:NumRawFiles; subplot(1,3,i); imagesc(PicCell{i}); axis([212 811 212 811]); end
title(['BaseLine = ' num2str(BaseLine)])
for PhantomIndex=1:1
[PicCell{PhantomIndex},Count(PhantomIndex),CountEw,LRFoutCell{PhantomIndex},SpectrumCell{PhantomIndex}]=...
    MexSPEngine_08insert_ANGER(LRF, pathFloodFile, LoopUC, PMTxy, SpectrumMax,...
    SpectrumWindow, LUT{PhantomIndex}, NumXlines, NumYlines, NumPM,PE,BaseLine);
disp(['Count: ',num2str(Count(PhantomIndex))]);
end
subplot(1,3,3); imagesc(PicCell{1}); axis([212 811 212 811]);
figname='Anger_XYF';
saveas(FigH,[WorkDir,'/Corr/Fig/PNG_H',num2str(Head,'%02i'),'_',figname,'.png'])
saveas(FigH,[WorkDir,'/Corr/Fig/FIG_H',num2str(Head,'%02i'),'_',figname,'.fig'])


toc; %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC
disp(' ');
disp('FIRST PART: Playing acquisitions with Statistic Method for spectrum check.');

% Play the acquisition %%%%%%%%%%%%%%%%%%%%%
SpectrumWindow=-1; % Which means that there is no energy window and LUTs are not used.
ShortLoop = 5e5;
PhantomIndex=1;

[PicCell{PhantomIndex},Count(PhantomIndex),CountEw,LRFoutCell{PhantomIndex},SpectrumCell{PhantomIndex}]=...
    MexSPEngine_08insert(LRF, path{PhantomIndex}, Loop(PhantomIndex), PMTxy, SpectrumMax,...
    SpectrumWindow, LUT{PhantomIndex}, NumXlines, NumYlines, NumPM,PE,BaseLine);
disp(['Count: ',num2str(Count(PhantomIndex))]);


% Figure for spectrum check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vector = 1:500;
Vector(:) = SpectrumCell{1}(1,1,:);
[tmp, PE] = max(Vector); PE = PE * 100;
SpectrumMax = ones(NumXlines,NumYlines) * PE;

figname='Spectrum';
FigH=figure(); set(FigH,'Position',[0 0 1920 1080]) ;
set(FigH,'NumberTitle','off','Name',figname);
plot(Vector)


disp('FIRST PART: Playing acquisitions with Statistic Method for line search method.');
SpectrumWindow=0.15; % Spectrum window is 2 * 15% = 30%
for PhantomIndex=1:NumRawFiles
    [PicCell{PhantomIndex},Count(PhantomIndex),CountEw,LRFoutCell{PhantomIndex},SpectrumCell{PhantomIndex}]=...
        MexSPEngine_08insert(LRF, path{PhantomIndex}, Loop(PhantomIndex), PMTxy, SpectrumMax,...
        SpectrumWindow, LUT{PhantomIndex}, NumXlines, NumYlines, NumPM,PE,BaseLine);
    disp(['Count: ',num2str(Count(PhantomIndex))]);
end

% % Figure for checking the acquisition play
% figure();
% for i=1:NumRawFiles subplot(NumRawFiles/2,2,i); imagesc(PicCell{i}); end
% return

toc; %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC
disp(' ');
disp('FIRST PART: Line search.');
disp(' ');


MidLine=MidLine1;
Extend='None'; % No restrictions on neither side.
Left=50;
Right=980;
% FindOnly='TRUE'; % 'FALSE' 'TRUE'
Tight='Both'; % No tight search on neither side.

saveas(FigH,[WorkDir,'/Corr/Fig/PNG_H',num2str(Head,'%02i'),'_',figname,'.png'])
saveas(FigH,[WorkDir,'/Corr/Fig/FIG_H',num2str(Head,'%02i'),'_',figname,'.fig'])
figname='Line search - default LRFs';
FigH=figure(); set(FigH,'Position',[0 0 1920 1080]); 
set(FigH,'NumberTitle','off','Name',figname);
colormap('pink');

NumLinesX1=0;NumLinesY1=0;NumLinesX2=0;NumLinesY2=0;

subplot(NumRawFiles/2,2,1);
NumLines=NumXlines; if (NumRawFiles==4) NumLines=NumX1Lines; end
while (NumLinesX1~=NumLines)
[LinesAllX1,NumLinesX1]=SnakesToRibbons_SP_C(PicCell{1},'X',MidLine,...
        Extend,Left,Right,'TRUE',Tight);
    if (NumLinesX1~=NumLines) MidLine=MidLine-1; end
end
[LinesAllX1,NumLinesX1]=SnakesToRibbons_SP_C(PicCell{1},'X',MidLine,...
        Extend,Left,Right,'FALSE',Tight);
%%

comdat('set','LineData', LinesAllX1);
comdat('set','ImageData', PicCell{1});

mo = corr_snakes();
waitfor(mo);

clear mo;
disp('continue');
LinesAllX1 = comdat('get','LineData');

%load('CorrectXSnk.mat')

%%
% !!!!!!!!!!!!!!!!!!!!!!!!
MidLine=MidLine2;
subplot(NumRawFiles/2,2,2);
NumLines=NumYlines; if (NumRawFiles==4) NumLines=NumY1Lines; end
while (NumLinesY1~=NumLines)
[LinesAllY1,NumLinesY1]=SnakesToRibbons_SP_C(PicCell{2},'Y',MidLine,...
        Extend,Left,Right,'TRUE',Tight);
    if (NumLinesY1~=NumLines) MidLine=MidLine-1; end  
end


[LinesAllY1,NumLinesY1]=SnakesToRibbons_SP_C(PicCell{2},'Y',MidLine,...
        Extend,Left,Right,'FALSE',Tight);

if (NumRawFiles==4)
    subplot(2,2,3);
    NumLines=NumX2Lines;
    while (NumLinesX2~=NumLines)
        [LinesAllX2,NumLinesX2]=SnakesToRibbons_SP_C(PicCell{3},'X',MidLine,...
            Extend,Left,Right,'TRUE',Tight);
        if (NumLinesX2~=NumLines) MidLine=MidLine-1; end
    end
    [LinesAllX2,NumLinesX2]=SnakesToRibbons_SP_C(PicCell{3},'X',MidLine,...
        Extend,Left,Right,'FALSE',Tight);
    
    subplot(2,2,4);
    NumLines=NumY2Lines;
    while (NumLinesY2~=NumLines)
        [LinesAllY2,NumLinesY2]=SnakesToRibbons_SP_C(PicCell{4},'Y',MidLine,...
            Extend,Left,Right,'TRUE',Tight);
        if (NumLinesY2~=NumLines) MidLine=MidLine-1; end
    end
    [LinesAllY2,NumLinesY2]=SnakesToRibbons_SP_C(PicCell{4},'Y',MidLine,...
        Extend,Left,Right,'FALSE',Tight);
end

%%
comdat('set','LineData', LinesAllY1);
comdat('set','ImageData', PicCell{2});

mo = corr_snakes();
waitfor(mo);

clear mo;
disp('continue');
LinesAllY1 = comdat('get','LineData');
%load('CorrectYSnk.mat')

%%

toc; %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC
disp(' ');
disp('FIRST PART: Generation of LUT table.');

if (NumRawFiles==4)
    [ShiftXY,Out]=A20_SplineToSnakes_08(PicCell{1},PicCell{2},LinesAllX1,LinesAllY1,LinesAllX2,LinesAllY2);
else
    [ShiftXY,Out]=A20_SplineToSnakes_08(PicCell{1},PicCell{2},LinesAllX1,LinesAllY1);
end

LUT=LUT_Generator_02(ShiftXY,LUT_mask1,NumRawFiles);


% SECOND PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plying the acquisitions with Statistic Method using LUT table.
% Spectrums to be generated in all line crosses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc; %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC
disp(' ');
disp('SECOND PART: Playing acquisitions with Statistic Method for generating SpectrumMax matrix.');

% Palying the acquisitions. %%%%%%%%%%%%%%%%%%%%%
SpectrumWindow=0.10; % Energy window = 2 * 10% = 20%

for PhantomIndex=1:NumRawFiles
    [PicCell{PhantomIndex},Count(PhantomIndex),CountEw,LRFoutCell{PhantomIndex},SpectrumCell{PhantomIndex}]=...
        MexSPEngine_08insert(LRF, path{PhantomIndex}, Loop(PhantomIndex), PMTxy, SpectrumMax,...
        SpectrumWindow, LUT{PhantomIndex}, NumXlines, NumYlines, NumPM,PE,BaseLine);
    disp(['Count: ',num2str(Count(PhantomIndex)),' Count in EnergyWindow: ',num2str(CountEw)]);
end


% % Figure for checking the acquisitions and the LUT table
% figure();
% for i=1:NumRawFiles subplot(NumRawFiles/2,2,i); imagesc(PicCell{i}); end


% Specturm check for all line crosses %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FigH=figure();
% colormap('pink');
% set(FigH,'NumberTitle','off','Name','Spectra in crosses for X acquisition');
% 
% for i=1:NumXlines
%     for j=1:NumYlines
%         Vector(:)=SpectrumCell{1}(i,j,:);
%         subplot(NumXlines,NumYlines,(j-1)*NumXlines+i);
%         plot(Vector);
%     end
% 
% end
% 
% 
% FigH=figure();
% colormap('pink');
% set(FigH,'NumberTitle','off','Name','Spectra in crosses for Y acquisition');
% 
% for i=1:NumXlines
%     for j=1:NumYlines
%         Vector(:)=SpectrumCell{2}(i,j,:);
%         subplot(NumXlines,NumYlines,(j-1)*NumXlines+i);
%         plot(Vector);
%     end
% 
% end


toc; %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC
disp(' ');
disp('SECOND PART: Running SpectrumMaxGenerator.');
SpectrumMax=SpectrumMaxGenerator(SpectrumCell,NumXlines,NumYlines,NumRawFiles); % Figyelni kell, hogy az illesztés ki van-e kommentelve a függvényben.


saveas(FigH,[WorkDir,'/Corr/Fig/PNG_H',num2str(Head,'%02i'),'_',figname,'.png'])
saveas(FigH,[WorkDir,'/Corr/Fig/FIG_H',num2str(Head,'%02i'),'_',figname,'.fig'])
figname='Spectrum max values in the crosses';
FigH=figure(); set(FigH,'Position',[0 0 1920 1080]); 
set(FigH,'NumberTitle','off','Name',figname);
colormap('pink');
imagesc(SpectrumMax);
disp(' ');

PE=mean(mean(SpectrumMax));


% THIRD PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Playing acquisitions with Statistic Method with using proper energy
% correction and windowing. The final LRF set is then received by smoothing
% the measured one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc; %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC
disp(' ');
disp('THIRD PART: Playing acquisitions with Statistic Method using energy correction and windowing for generating final LRF set.');


% Running play acquisition. %%%%%%%%%%%%%%%%%%%%%
SpectrumWindow=0.10; % Energy window: 2 * 10% = 20%.

for PhantomIndex=1:NumRawFiles
    [PicCell{PhantomIndex},Count(PhantomIndex),CountEw,LRFoutCell{PhantomIndex},SpectrumCell{PhantomIndex}]=...
        MexSPEngine_08insert(LRF, path{PhantomIndex}, Loop(PhantomIndex), PMTxy, SpectrumMax,...
        SpectrumWindow, LUT{PhantomIndex}, NumXlines, NumYlines, NumPM,PE,BaseLine);
    disp(['Count: ',num2str(Count(PhantomIndex)),' Count in EnergyWindow: ',num2str(CountEw)]);
end

% figure();
% for i=1:NumRawFiles subplot(NumRawFiles/2,2,i); imagesc(PicCell{i}); end


toc; %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC
disp(' ');
% disp('THIRD PART: Smoothing LRF set: Generating final LRF set.');
disp('THIRD PART - ORIGINALLY: Smoothing LRF set: Generating final LRF set, but not for Clinical case.');
% disp('no LRF calc')

% Final smoothed sqrt(LRF)s are generated. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Filter1=0.05;
% 
% [LRF,LRFoutS]=SmoothLRF_07_par(LRFoutCell,0,NumXlines,NumYlines,NumPM,NumRawFiles, Filter1,Dx,dx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


toc; %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC

% LRFstartsum=reshape(sum(LRFstart,1),1024,1024);
% figure
% subplot(1,2,1)
% imagesc(LRFstartsum)
% title('start lrf')
% LRFsum=reshape(abs(sum(LRF,1)),1024,1024);
% subplot(1,2,2)
% imagesc(LRFsum)
% title('act lrf')

% FOURTH PART %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating all correction tables.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(' ');
disp('FOURTH PART: Playing acquisitions with Statistic Method using final LRFs. - In Clinical case using starting LRFs.');

% Playing acquisitions %%%%%%%%%%%%%%%%%%%%%
SpectrumWindow=0.10; % Energy window: 2 * 10% = 20%

for PhantomIndex=1:NumRawFiles
    [PicCell{PhantomIndex},Count(PhantomIndex),CountEw,LRFoutCell{PhantomIndex},SpectrumCell{PhantomIndex}]=...
        MexSPEngine_08insert(LRF, path{PhantomIndex}, Loop(PhantomIndex), PMTxy, SpectrumMax,...
        SpectrumWindow, LUT{PhantomIndex}, NumXlines, NumYlines, NumPM,PE,BaseLine);
    disp(['Count: ',num2str(Count(PhantomIndex))]);
end

% figure();
% colormap('pink');
% for i=1:NumRawFiles subplot(NumRawFiles/2,2,i); imagesc(-PicCell{i},[0,-min(min(PicCell{i}(500:700,300:700)))]); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Removing too light side pixel rows and columns.
for i=1:NumRawFiles;
    Vector=1:1024;
    NonZeroI = ( sum(PicCell{i},2)~=0 );
    NonZeroI = NonZeroI .* Vector';
    NonZeroItrunc=NonZeroI(NonZeroI>0);
    minI=min(NonZeroItrunc); maxI=max(NonZeroItrunc);
    PicCell{i}(minI,:)=0; PicCell{i}(maxI,:)=0;
    
    
    NonZeroJ = ( sum(PicCell{i},1)~=0 );
    NonZeroJ = NonZeroJ .* Vector;
    NonZeroJtrunc=NonZeroJ(NonZeroJ>0);
    minJ=min(NonZeroJtrunc); maxJ=max(NonZeroJtrunc);
    PicCell{i}(:,minJ)=0; PicCell{i}(:,maxJ)=0;
end



toc; %TOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOCTOC

% Line search.
disp(' ');
disp('FOURTH PART: Line search.');

MidLine=MidLine1; % A középvonal mentén keresünk
Extend='None'; % Nem kényszerítünk pontokat a jobb és bal szélekre.
Left=50;
Right=980;
FindOnly='FALSE'; % 'FALSE' 'TRUE'
Tight='Both'; % Egyik oldalon sem alkalmazunk tight searchot.

saveas(FigH,[WorkDir,'/Corr/Fig/PNG_H',num2str(Head,'%02i'),'_',figname,'.png'])
saveas(FigH,[WorkDir,'/Corr/Fig/FIG_H',num2str(Head,'%02i'),'_',figname,'.fig'])
figname='Line search - actual LRFs';
FigH=figure(); set(FigH,'Position',[0 0 1920 1080]); 
set(FigH,'NumberTitle','off','Name',figname);
colormap('pink');

subplot(NumRawFiles/2,2,1);
[LinesAllX1,NumLinesX1]=SnakesToRibbons_SP_C(PicCell{1},'X',MidLine,...
        Extend,Left,Right,FindOnly,Tight);
%%

comdat('set','LineData', LinesAllX1);
comdat('set','ImageData', PicCell{1});

mo = corr_snakes();
waitfor(mo);

clear mo;
disp('continue');
LinesAllX1 = comdat('get','LineData');

%%

MidLine=MidLine2; % A középvonal mentén keresünk
subplot(NumRawFiles/2,2,2);
[LinesAllY1,NumLinesY1]=SnakesToRibbons_SP_C(PicCell{2},'Y',MidLine,...
        Extend,Left,Right,FindOnly,Tight);
%%

comdat('set','LineData', LinesAllY1);
comdat('set','ImageData', PicCell{2});

mo = corr_snakes();
waitfor(mo);

clear mo;
disp('continue');
LinesAllY1 = comdat('get','LineData');

%%
if (NumRawFiles==4)
    subplot(2,2,3);
    
    [LinesAllX2,NumLinesX2]=SnakesToRibbons_SP_C(PicCell{3},'X',MidLine,...
        Extend,Left,Right,FindOnly,Tight);
    subplot(2,2,4);
    [LinesAllY2,NumLinesY2]=SnakesToRibbons_SP_C(PicCell{4},'Y',MidLine,...
        Extend,Left,Right,FindOnly,Tight);
end

% Generating linearity tables.
disp(' ');
disp('FOURTH PART: Generating linearity correction tables.');


[ShiftXY,Out]=A20_SplineToSnakes_08(PicCell{1},PicCell{2},LinesAllX1,LinesAllY1);

NumIlines=size(ShiftXY,1);
NumJlines=size(ShiftXY,2);


HalfIline=(NumIlines+1)/2; HalfJline=(NumJlines+1)/2;
[ThetaX,ThetaY,Out]=SplineToCrosses(ShiftXY,HalfIline,HalfJline, NumIlines, NumJlines,Dx,dx);

[LinX,LinY]=LoadXYtr_Mex(ThetaX,ThetaY,HalfIline,HalfJline, NumIlines, NumJlines,Dx,dx);

toc;

% Generating energy correction table.
disp(' ');
disp('FOURTH PART: Generating energy correction table.');

ECcore=imresize(SpectrumMax,Dx/dx/2); ECcore=imresize(ECcore,2);

EC = zeros(1024,1024);

[size1,size2]=size(ECcore);

EC( (512-size1/2+1):(512+size1/2),(512-size2/2+1):(512+size2/2) ) = ECcore(:,:);

EC=EC';

saveas(FigH,[WorkDir,'/Corr/Fig/PNG_H',num2str(Head,'%02i'),'_',figname,'.png'])
saveas(FigH,[WorkDir,'/Corr/Fig/FIG_H',num2str(Head,'%02i'),'_',figname,'.fig'])
figname='Energy table';
FigH=figure(); set(FigH,'Position',[0 0 1920 1080]); 
set(FigH,'NumberTitle','off','Name',figname);
colormap('pink');
imagesc(EC,max(max(EC))*[0.85,1.0]);


mask=zeros(1024,1024);

maskI=round((NumYlines-1)*Dx/dx/2); maskJ=round((NumXlines-1)*Dx/dx/2);
mask( (512-maskI+1):(512+maskI),(512-maskJ+1):(512+maskJ) ) = 1;



% for i=1:length(pathFloodFile)
%     if ( pathFloodFile(i)=='\' )
%         pathFloodFile(i)='/';
%     end;
% end



% Playing flood acquisition.
disp(' ');
disp('FOURTH PART: Playing flood acquisition with Statistic Method using linearity and energy correction for generating uniformity calibration.');

SpectrumWindow=0.10;

UC=ones(1024);

[UC,Count,CountEw]=...
    MexSPEngine_10insertUCECLin(LRF, pathFloodFile, LoopUC, PMTxy,...
    NumPM, SpectrumWindow,EC,UC,LinX,LinY,PE,BaseLine);
disp(['Count: ',num2str(Count),' Count in EnergyWindow: ',num2str(CountEw)]);


% Generating uniformity correction table.
disp(' ');
disp('FOURTH PART: Generating uniformity correction table.');

UC=imresize(UC,0.25,'bilinear');
UC=imresize(UC,2,'bilinear'); UC=imresize(UC,2,'bilinear');
UC(UC==0)=1;
UC=ones(1024,1024)./UC;
UC=UC.*mask;

UC = UC / sum(sum(UC)) * sum(sum(UC>0));
saveas(FigH,[WorkDir,'/Corr/Fig/PNG_H',num2str(Head,'%02i'),'_',figname,'.png'])
saveas(FigH,[WorkDir,'/Corr/Fig/FIG_H',num2str(Head,'%02i'),'_',figname,'.fig'])
figname='Uniformity table';
FigH=figure(); set(FigH,'Position',[0 0 1920 1080]); 
set(FigH,'NumberTitle','off','Name',figname);
colormap('pink');
imagesc(UC);
caxis([0 0.01])



% FINISHING... %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replaying acquisitions with all corrections.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Re-playing acquisitions.
disp(' ');
disp('FINISHING: Re-playing acquisitions.');

% UC=ones(1024,1024);
SpectrumWindow=0.10;
for PhantomIndex=1:NumRawFiles
    [PicCell{PhantomIndex},Count(PhantomIndex),CountEw]=...
        MexSPEngine_10insertUCECLin(LRF, path{PhantomIndex}, Loop(PhantomIndex), PMTxy,...
        NumPM, SpectrumWindow,EC,UC,LinX,LinY,PE,BaseLine);
    disp(['Count: ',num2str(Count(PhantomIndex)),' Count in EnergyWindow: ',num2str(CountEw)]);
end

saveas(FigH,[WorkDir,'/Corr/Fig/PNG_H',num2str(Head,'%02i'),'_',figname,'.png'])
saveas(FigH,[WorkDir,'/Corr/Fig/FIG_H',num2str(Head,'%02i'),'_',figname,'.fig'])
figname='X and Y acquisitions - replay';
FigH=figure(); set(FigH,'Position',[0 0 1920 1080]); 
set(FigH,'NumberTitle','off','Name',figname);
colormap('pink');
for i=1:NumRawFiles subplot(NumRawFiles/2,2,i); imagesc(PicCell{i}); end

SpectrumWindow=0.10;
[Pic,Count,CountEw]=...
    MexSPEngine_10insertUCECLin(LRF, pathFloodFile, LoopUC, PMTxy,...
    NumPM, SpectrumWindow,EC,UC,LinX,LinY,PE,BaseLine);
disp(['Count: ',num2str(Count),' Count in EnergyWindow: ',num2str(CountEw)]);

saveas(FigH,[WorkDir,'/Corr/Fig/PNG_H',num2str(Head,'%02i'),'_',figname,'.png'])
saveas(FigH,[WorkDir,'/Corr/Fig/FIG_H',num2str(Head,'%02i'),'_',figname,'.fig'])
figname='Flood acquisition - replay';
FigH=figure(); set(FigH,'Position',[0 0 1920 1080]); 
set(FigH,'NumberTitle','off','Name',figname);
colormap('pink');
imagesc(Pic);


% Saving corrections.
disp(' ');
disp('FINISHING: Saving corrections.');

if ~isdir([WorkDir,'/Corr'])
    mkdir([WorkDir,'/Corr']);
end

saveas(FigH,[WorkDir,'/Corr/Fig/PNG_H',num2str(Head,'%02i'),'_',figname,'.png'])
saveas(FigH,[WorkDir,'/Corr/Fig/FIG_H',num2str(Head,'%02i'),'_',figname,'.fig'])

save([WorkDir,'/Corr/CorrHead',num2str(Head,'%02i'),'.mat'],...
     'LRF','LinX','LinY','UC','EC','PMTxy','PE','BaseLine');
close all
toc;
end
