function Fnc_RefreshUnifCorr(Head, NumXlines, NumYlines, Dx, dx, ...
             NumPM, pathXfile, pathYfile, pathFloodFile, WorkDir)


tic; %TICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTICTIC

% Loading correction data file
load([WorkDir,'/Corr/CorrHead',num2str(Head,'%02i'),'.mat'],...
     'LRF','LinX','LinY','UC','EC','PMTxy', 'BaseLine', 'PE');

% Getting length of flood file (LoopUC)
StreamFile=fopen(pathFloodFile,'r');
fseek(StreamFile,0,'eof');
LoopUC=floor(ftell(StreamFile)/(72*2+4*2));
fclose(StreamFile);
% Value is rounded by 10,000
LoopUC=LoopUC-rem(LoopUC,1e4);

% Backslash <-> slash; Path prepared for C function
pathFloodFile = regexprep(pathFloodFile,'\','/');


% Initialization...
SpectrumWindow=0.10;

UC=ones(1024);

PE = mean( EC( EC~=0 ) );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Playing flood acquisition
disp('Playing flood acquisition...');
[UC,Count,CountEw]=...
    MexSPEngine_10insertUCECLin(LRF, pathFloodFile, LoopUC, PMTxy,...
    NumPM, SpectrumWindow,EC,UC,LinX,LinY,PE,BaseLine);
disp(['Count: ',num2str(Count),' CountEw: ',num2str(CountEw)]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formatting UC table

% Squeezing and pulling the UC image by a factor of 4.
UC=imresize(UC,0.25,'bilinear');
UC=imresize(UC,2,'bilinear'); UC=imresize(UC,2,'bilinear');

% Generating reciprocal
UC(UC==0)=1;
UC=ones(1024,1024)./UC;

% Creating a mask in the rectangle determined by the phantom lines of
% linearity phantom.
mask=zeros(1024,1024);
maskI=round((NumYlines-1)*Dx/dx/2); maskJ=round((NumXlines-1)*Dx/dx/2);
mask( (512-maskI+1):(512+maskI),(512-maskJ+1):(512+maskJ) ) = 1;

% Masking the UC table
UC=UC.*mask;

% Normalizing UC table
UC = UC / sum(sum(UC)) * sum(sum(UC>0));
FigH=figure();
colormap('pink');
set(FigH,'NumberTitle','off','Name','Uniformity table');
imagesc(UC);
caxis([0 0.01])
disp('Uniformity correction table created.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Playing back flood acquisition for testing UC table
disp('Playing back flood acquisition...');
[Pic,Count,CountEw]=...
    MexSPEngine_10insertUCECLin(LRF, pathFloodFile, LoopUC, PMTxy,...
    NumPM, SpectrumWindow,EC,UC,LinX,LinY,PE,BaseLine);
disp(['Count: ',num2str(Count),' CountEw: ',num2str(CountEw)]);

FigH=figure();
set(FigH,'NumberTitle','off','Name','Flood acquisition - replay');
colormap('pink');
imagesc(Pic);


% Dialog: save or not save?
button = questdlg('Would you like to save the new UC table?','Confirmation','Yes','No','Yes');

if strcmp(button,'Yes')
    save([WorkDir,'/Corr/CorrHead',num2str(Head,'%02i'),'.mat'],...
         'LRF','LinX','LinY','UC','EC','PMTxy');
    disp('New uniformity correction table saved.');
else
    disp('New uniformity correction table table has not been saved!');
end

EllapsedTime = toc;
 
disp([ 'Ellapsed time: ', num2str(EllapsedTime) ]);
