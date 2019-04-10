function Fnc_PlayAcquisition(handles)

% //-> Itt van bevarrva a fejek szama.

tic

NumHeads = handles.NumHeads;

if NumHeads ~= 10
    disp('Play acquisition has not been generalized yet. Number of heads must be 10');
    disp('new: return was commented')
%     return;
end


WorkDir = handles.WorkDir;
targetDir = handles.targetDir;

GeneralLoop = 0; %1e5;

for h=1:NumHeads
    
    pathRaw = [targetDir,'\digitalraw_H',num2str(h,'%02i'),'.dat'];
    pathCorr = [WorkDir, '\Corr\CorrHead',num2str(h,'%02i'),'.mat'];
    load(pathCorr) %, 'EC', 'LinX', 'LinY', 'LRF', 'PMTxy', 'UC' , 'BaseLine', 'PE');
        
      
    % Loop beállítása
    StreamFile=fopen(pathRaw,'r');
    fseek(StreamFile,0,'eof');
    Loop=floor(ftell(StreamFile)/(72*2+4*2));
    fclose(StreamFile);

    if ( GeneralLoop~=0 && GeneralLoop<Loop )
        Loop=GeneralLoop;
    end
    
    Loop=Loop-rem(Loop,1e4);
    
    SpectrumWindow=0.1;
%     PE=17000;
 
    [PicCell{h},Count(h),CountEw(h),SPEC1{h},SPEC2{h}]=...
        MexSPEngine_10insertUCECLin( LRF, pathRaw, Loop,  PMTxy,...
        handles.NumPM, SpectrumWindow, EC, UC, LinX, LinY, PE, BaseLine); %NEW: SPEC1 and SPEC2
    disp(['Count: ',num2str(Count(h)),' CountEw: ',num2str(CountEw(h))]);
     
    figure();  
    colormap('hot');    
    imagesc(PicCell{h})
    
    %NEW
    figure();
    plot(1:numel(SPEC1{h}),smooth(SPEC1{h},40),'linewidth',1.5) 
    hold on
    plot(1:numel(SPEC2{h}),smooth(SPEC2{h},40),'linewidth',1.5) 
    legend('spectrum','photopeak selection')
    title(strcat('Node',32,num2str(h)))
    disp(strcat('counts:',32,num2str(sum(SPEC1{h})),32,'- photopeak counts:',32,num2str(sum(SPEC2{h})),32,'- node:',32,num2str(h)))
    %NEW
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading template DCM file
% try
%     dcmFile = 'DCMnanoTemplate.dcm';
%     newDcmFile=[targetDir,'\new.dcm'];
%     
%     header = dicominfo( dcmFile );
%     
%     image = dicomread( dcmFile );
% catch ME
%    disp(ME.message); 
% end
% 
% if size(image) ~= [160 160 1 10]
%     disp('Play acquisition has not been generalized yet. Template DICOM dimensions must be [160 160 1 10]');
%     return;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Formatting Pics
dx=handles.dx;

imagesize = round((0.2/dx) * 480 /2 ) * 2 ;
image = zeros(imagesize,imagesize,10);
% image_bin2x2 = zeros(240,240,1,10);
slice = zeros(imagesize,imagesize);

for h=1:NumHeads %NumHeads
    
    slice(:,:)=PicCell{h}( (512-imagesize/2+1 : 512+imagesize/2),(512-imagesize/2+1 : 512+imagesize/2) );
    slice = flip( slice, 2 );
%     slice160 = (imresize(slice,[160,160],'box')*9); % 9: mert 3X3-anként összevontuk a pixeleket.
%     slice240 = (imresize(slice,[240,240],'box')*4); % 9: mert 3X3-anként összevontuk a pixeleket.
    image(:,:,h) = slice(:,:);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create DICOM with inline Python call
TotCount = {Count, CountEw};
save([handles.targetDir,'\CountsAndEWC.mat'],'TotCount');

save([handles.targetDir,'\corrected_images.mat'],'image');

% cmdstr = ['START /WAIT ', 'python _PyDicom\CreateDICOMfromMat.py ', handles.targetDir,'\dicom.mat ', '_TMP/DCMnanoTemplate.dcm '...
%     , handles.targetDir ,'\PythonDCM.dcm' ]

% cmdstr = ['START /WAIT ', 'python _PyDicom\CreateDICOMfromMat.py ', handles.targetDir,'\dicom.mat ', '_TMP/DCMnanoTemplate.dcm '...
%     , handles.targetDir ,'\PythonDCM.dcm' ]
% s = system(cmdstr);

toc;

end