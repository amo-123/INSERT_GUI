function image = Fnc_PlayAcquisition_n(CorrDir,digitalrawDir,saveimageFilePath,n_ang)

% digitalrawDir = '';
% saveimageFilePath = [digitalrawDir,'\corrected_images.mat'];

%n_ang = 28;                                                                % number of angles

p1 = [1.19291,-0.1237386];
p2 = [0.717933,8.454];

tic

NumHeads = 10;
NumPM=72;
dx=0.2;
SpectrumWindow=0.1;
GeneralLoop = 0; 
PicCell=cell(1,10);
imsize = round((0.2/dx) * 480 /2 ) * 2 ;

image = zeros(imsize,imsize,n_ang,NumHeads);
slice = zeros(imsize,imsize);

aa = zeros(imsize,imsize,n_ang);

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
    
    Loop0 = Loop;
    for i_ang=1:n_ang
        
        if ( i_ang <= 18 )
            ff = ( p1(1) * i_ang + p1(2) ) / n_ang;
        else
            ff = ( p2(1) * i_ang + p2(2) ) / n_ang;
        end
        ff = max([0,ff]);  ff = min([ff,1]);
        
        Loop = floor( Loop0 * ff );
        Loop=Loop-rem(Loop,1e4);
    
        [PicCell{h},Count,CountEw]=...
            MexSPEngine_10insertUCECLin( LRF, pathRaw, Loop,  PMTxy,...
            NumPM, SpectrumWindow, EC, UC, LinX, LinY, PE, BaseLine);
        disp(['Count: ',num2str(Count),' CountEw: ',num2str(CountEw)]);
        
        slice(:,:)=PicCell{h}( (512-imsize/2+1 : 512+imsize/2),(512-imsize/2+1 : 512+imsize/2) );
        slice = flip( slice, 2 ); 
        image(:,:,i_ang,h) = slice(:,:);
        
    end
    
    figure();
    
    colormap('pink');
    
    imagesc(PicCell{h})

end

for h=1:NumHeads
    for i_ang=n_ang:-1:2
        image(:,:,i_ang,h) = image(:,:,i_ang,h) - image(:,:,i_ang-1,h);
    end
end

save(saveimageFilePath,'image'); 

toc;

end