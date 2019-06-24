%% Extract Single Node
%Code to extract and visualize the Frame of one Node only, run Main_Split first to get FRAME_NODES variable

clearvars -except filename filepath modality FRAME_NODE num_nodes
clc

%% Inzializations

flag.zeros_saturations_filtering = 0;
flag.offsets_subtraction = 0;

valori_impulsazione=100:100:900;
zero_thr=10;
saturation_thr=4085;
en_window_perc_width=0.05; 

%% Geometry

if strcmp(modality,'Preclinical')
    load('INSERT_Preclinical.mat');
elseif strcmp(modality,'Clinical')
    load('INSERT_Clinical.mat');
end

%% Node selection

NODE=menu('Select the Node:',num2cell(reshape(1:num_nodes,1,num_nodes)));
Frame=FRAME_NODE{NODE,1};
disp(strcat('Node:',32,num2str(NODE)))

%% Visualization

%channels histograms
disp('Single channels histograms')
histograms=PlotHistograms(modality,Frame,NODE,zero_thr);

%spectrum
bins=1000;
[spectrum,spectrum_filt]=PlotSpectrum(Frame,NODE,bins,saturation_thr);

%cog image
%-peak selection
fitting=monoG(spectrum_filt.x,spectrum_filt.y,'Select photopeak position',0.8,0.8);
x_peak=round(fitting.b1); %mean
FWHM=fitting.c1/sqrt(2)*2.355;
RES_nocalib=FWHM/x_peak;
%-energy windowing
Filt.E_min=round((1-en_window_perc_width).*x_peak);
Filt.E_max=round((1+en_window_perc_width).*x_peak);
vline_spectrum_plot(spectrum_filt.y,Filt.E_min,Filt.E_max)
%-baseline subtraction for mcog
Filt.baseline=input(strcat('Baseline value =',32)); %[ADCchannels] - 600(normal)|500(after offset subtraction)|200
figure('units','normalized','outerposition',[0 0 1 1])
%-mcog rec
[IMG_cog.x_rec,IMG_cog.y_rec,IMG_cog.Counts,IMG_cog.energy_window,IMG_cog.number_det_ch]=CentroidReconstruction(Frame,Par,Filt);

%% Elaboration

if ~exist('Frame_original','var')
    Frame_original=Frame;
end

%zeros and saturations filtering
if flag.zeros_saturations_filtering == 1   
    f = figure('units','normalized','outerposition',[0 0 1 1],'Name',filename,'NumberTitle','off');
    [~]=filt_zeros_saturations(Frame,zero_thr,saturation_thr,'plot');
    choice=menu('Do you want to filter all the nodes?','Yes - only saturations','Yes - zeros and saturations', 'No');
    if choice==1
        [Frame]=filt_zeros_saturations(Frame,-1,saturation_thr,'filt');
        disp('all nodes filtered - only saturations')
    elseif choice==2
            [Frame]=filt_zeros_saturations(Frame,zero_thr,saturation_thr,'filt');
            disp('All nodes filtered - zeros and saturations')
    else disp ('No filtering')
    end
end

%offset subtraction
if flag.offsets_subtraction == 1
    [offsets,coeff]=fcn_calcolo_offsets_ch(valori_impulsazione);
    if isvector(offsets)
        if numel(offsets) ~= size(Frame,2)
            disp('ERROR! Vector size mismatch')
        else
            if iscolumn(offset)
                offset = offset';
            end
        end
        Frame = Frame + repmat(-offsets, size(Dataset,1),1);
        disp('Offsets subtraction performed')
    end
end










