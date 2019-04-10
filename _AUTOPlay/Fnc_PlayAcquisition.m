function image = Fnc_PlayAcquisition(CorrDir,digitalrawDir,saveimageFilePath)

% digitalrawDir = '';
% saveimageFilePath = [digitalrawDir,'\corrected_images.mat'];

tic

NumHeads = 10;
NumPM=72;
dx=0.2;
SpectrumWindow=0.1;
GeneralLoop = 0; 
PicCell=cell(1,10);

for h=1:NumHeads
    
    pathRaw = [digitalrawDir,'\digitalraw_H',num2str(h,'%02i'),'.dat'];
    pathCorr = [CorrDir '\CorrHead',num2str(h,'%02i'),'.mat'];
    load(pathCorr) %, 'EC', 'LinX', 'LinY', 'LRF', 'PMTxy', 'UC' , 'BaseLine', 'PE');

    StreamFile=fopen(pathRaw,'r');
    fseek(StreamFile,0,'eof');
    Loop=floor(ftell(StreamFile)/(72*2+4*2));
    fclose(StreamFile);

    if ( GeneralLoop~=0 && GeneralLoop<Loop )
        Loop=GeneralLoop;
    end
    
    Loop=Loop-rem(Loop,1e4);
    
    [PicCell{h},Count,CountEw]=...
        MexSPEngine_10insertUCECLin( LRF, pathRaw, Loop,  PMTxy,...
        NumPM, SpectrumWindow, EC, UC, LinX, LinY, PE, BaseLine);
    disp(['Count: ',num2str(Count),' CountEw: ',num2str(CountEw)]);
    
    
    figure();
    
    colormap('pink');
    
    imagesc(PicCell{h})
    
end


imagesize = round((0.2/dx) * 480 /2 ) * 2 ;
image = zeros(imagesize,imagesize,10);
slice = zeros(imagesize,imagesize);

for h=1:NumHeads 
    slice(:,:)=PicCell{h}( (512-imagesize/2+1 : 512+imagesize/2),(512-imagesize/2+1 : 512+imagesize/2) );
    slice = flip( slice, 2 );
    image(:,:,h) = slice(:,:);    
end

save(saveimageFilePath,'image'); 

toc;

end