function [OFFSETS,COEFF] = fcn_calcolo_offsets_ch(valori_impulsazione)

%load function
FilterSpec = '*.data';
[filename,filepath] = uigetfile([pwd,'\',FilterSpec], 'Select impulsation .data file', 'MultiSelect', 'off');
[Frame,Node,~,modality]=openDataFile(filename,filepath);

%load proper reorder array and geometrical reorder
[Phys_order]=INSERT_reorder(modality);
Frame(:,Phys_order) = Frame;
Node = round(Node);
%calculate number of nodes
num_nodes = max(Node);

for n=1:num_nodes 
    
    disp(strcat('Current Node:',32,num2str(n)))
    Frame_node=Frame(Node==n,:);
    [~,coefficienti]=new_equalizzazioneCanali(valori_impulsazione,Frame_node);
    
    [offsets] = fcn_Extract_Offsets(coefficienti);
    OFFSETS.(strcat('node',num2str(n)))=offsets;
    COEFF.(strcat('node',num2str(n)))=coefficienti;
    
end   
end

