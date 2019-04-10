
for imodule=1:20;
    
    tic
    pathlrfs=['C:\btolgyesi\POSTA\INSERT\clinical\Elaboration\module' num2str(imodule) '_Co57_gain15_th30_HV35e4_flood__LRF.mat' ];
    [LRF,PMTxy]=extraLRFwriteFromSfit(dx,NumPM,pathlrfs);
    
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
disp('FIRST PART: Playing acquisitions with Statistic Method for spectrum check.');

% Play the acquisition %%%%%%%%%%%%%%%%%%%%%
SpectrumWindow=-1; % Which means that there is no energy window and LUTs are not used.
ShortLoop = 5e5;
PhantomIndex=1;

[PicCell{PhantomIndex},Count(PhantomIndex),CountEw,LRFoutCell{PhantomIndex},SpectrumCell{PhantomIndex}]=...
    MexSPEngine_08insert(LRF, path{PhantomIndex}, min(Loop(PhantomIndex), ShortLoop ), PMTxy, SpectrumMax,...
    SpectrumWindow, LUT{PhantomIndex}, NumXlines, NumYlines, NumPM,PE);
disp(['Count: ',num2str(Count(PhantomIndex))]);


% Figure for spectrum check %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vector = 1:500;
Vector(:) = SpectrumCell{1}(1,1,:);
[tmp, PE] = max(Vector); PE = PE * 100;
SpectrumMax = ones(NumXlines,NumYlines) * PE;


FigH=figure(1588);
subplot(4,5,imodule)
set(FigH,'NumberTitle','off','Name','Spectrum');
plot(Vector)


disp('FIRST PART: Playing acquisitions with Statistic Method for line serach method.');
SpectrumWindow=0.15; % Spectrum window is 2 * 15% = 30%
for PhantomIndex=1:NumRawFiles
    [PicCell{PhantomIndex},Count(PhantomIndex),CountEw,LRFoutCell{PhantomIndex},SpectrumCell{PhantomIndex}]=...
        MexSPEngine_08insert(LRF, path{PhantomIndex}, Loop(PhantomIndex), PMTxy, SpectrumMax,...
        SpectrumWindow, LUT{PhantomIndex}, NumXlines, NumYlines, NumPM,PE);
    disp(['Count: ',num2str(Count(PhantomIndex))]);
end

% % Figure for checking the acquisition play
  h=figure(imodule);
 for i=1:NumRawFiles 
     subplot(1,3,i); 
     imagesc(PicCell{i}); 
 end
 subplot(1,3,3); 
     imagesc(reshape(sum(LRF,1),1024,1024));
     
      
     title(num2str(imodule))
     saveas(h,['C:\btolgyesi\POSTA\INSERT\insert_calibgui_modpz\_GenerateCorrections\figH05p\' num2str(imodule) '.fig'])
     


toc
end