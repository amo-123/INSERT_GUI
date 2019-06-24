function [ x_rec,y_rec,Counts,energy_window ,number_det_ch] = CentroidReconstruction( Frame,Par,Filt)

%Define number of pixels in X and Y for image sampling
pixel = 0.2;%mm
biny = ceil(Par.cryst_lung_y/pixel);
binx = ceil(Par.cryst_lung_x/pixel);
%% Definizione variabili
%Copy energy filter values
E_min=Filt.E_min;
E_max=Filt.E_max;
%Baseline definition
if isfield(Filt,'baseline')
    Baseline = Filt.baseline;
else
    Baseline = mean(mean(Frame)); %ADCchannels
end
%Detection channels centers
X_mat=repmat(Par.X_c,size(Frame,1),1);
Y_mat=repmat(Par.Y_c,size(Frame,1),1);

%% Centroid method
%Baseline subraction
Frame_sub=Frame-Baseline;
Frame_sub(Frame_sub<0)=0;

%Number of active detection channels per event, after baseline subtraction
number_det_ch=sum(Frame_sub>0,2);

%Centroid Formula
Energy_sub=sum(Frame_sub,2);
x_rec=sum((Frame_sub.*X_mat),2)./Energy_sub;
y_rec=sum((Frame_sub.*Y_mat),2)./Energy_sub;

%% Filters
%energy window
Energy = sum(Frame,2);
F1=Energy>E_min;
F2=Energy<E_max;
%Coordinates reconstructed in the Field of View (FOV)
F_tempX = ((x_rec>=0)&(x_rec<=Par.cryst_lung_x))>0;
F_tempY = ((y_rec>=0)&(y_rec<=Par.cryst_lung_y))>0;
F_FOV = F_tempX.*F_tempY;
%Filter dependent on the number of detection channels activated (influences only the image)
F3=number_det_ch<=size(Frame,2);
F4=number_det_ch>=3;

%Image filter (for Centroid reconstruction image)
filter_image=find((F1.*F2.*F3.*F4.*F_FOV)==1);

y_image=y_rec(filter_image,:)';
x_image=x_rec(filter_image,:)';
%two "fake" events are added at the corners of the FOV to ensure the histogram to get all the FOV (they are not considered in the output)
y_image  = [y_image 0 Par.cryst_lung_y];
x_image  = [x_image 0 Par.cryst_lung_x];
[Counts,Centers]=hist3([y_image' x_image'],[biny binx]);
% Counts(1,1)=Counts(1,1)-1;
% Counts(end,end)=Counts(end,end)-1;

%% Image
% immagine=flipud(immagine);
imagesc([0 Par.cryst_lung_x],[0 Par.cryst_lung_y],Counts)
hold on
set(gca,'Fontsize',15,'Fontname','TimesNewRoman','FontWeight','bold')
xlabel('x [mm]')
ylabel('y [mm]')
hold on
colormap hot
axis equal
axis tight
colorb=colorbar;
set(colorb,'Fontsize',12);
set(colorb,'FontWeight','bold');
set(colorb,'Fontname','Arial');
drawnow;

%% Output reconstruction
%Filter energy window (for LRF estimation)
energy_window=find((F1.*F2.*F_FOV)==1);
y_rec=y_rec(energy_window,:)';
x_rec=x_rec(energy_window,:)';

end

