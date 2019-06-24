function [output]=PlotImages(Frame,modality,Filt,num_nodes,Node,FName)

if strcmp(modality,'Preclinical')
    load('INSERT_Preclinical.mat');
elseif strcmp(modality,'Clinical')
    load('INSERT_Clinical.mat');
end

clc
disp('Centroid reconstruction')
% figure('units','normalized','outerposition',[0 0 1 1])

figure('Name',FName,'NumberTitle','off','units','normalized','outerposition',[0 0 1 1]);
cols = ceil(sqrt(num_nodes));
rows= ceil(num_nodes/cols);
for n = 1:num_nodes
    Filt_node = Filt;
    Frame_node = Frame(Node == n,:);
    if iscell(Filt.E_min)
        Filt_node.E_min=Filt.E_min{n};
        Filt_node.E_max=Filt.E_max{n};
    end
    %Filt_node.baseline=Filt.baseline{n};
    subplot(rows,cols,n)
    [ x_rec, y_rec, Counts, energy_window, number_det_ch ] = CentroidReconstruction( Frame_node,Par,Filt_node);
    title(['Node ',num2str(n)])
    
    output.(['node_',num2str(n)]).x_rec = x_rec;
    output.(['node_',num2str(n)]).y_rec = y_rec;
    output.(['node_',num2str(n)]).Image = Counts;
    output.(['node_',num2str(n)]).En_filter = energy_window; % boolean vector (high values defines the events that survived filtering)
    output.(['node_',num2str(n)]).Active_Channels = number_det_ch; %number of active channels after baseline subtraction
end


end