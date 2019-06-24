function [Data,nodeID,timestamp,modality] = openDataFile(File_Name,File_Path,num_events)

disp(strcat('Loading .data File:',{' '},File_Name));

%initializations and parameters
time_res = 25e-9; 
disp(strcat('time resolution =',32,num2str(time_res)));

%% File opening

fid = fopen([File_Path,File_Name],'r','b'); %windows

%check file dimension
fseek(fid,0,'eof');
FileSize=ftell(fid); %location of the file position indicator is indicated in bytes from the beginning of the file
fseek(fid,0,'bof');

file = fread(fid,8,'uint16=>uint16','b');
fseek(fid,0,'bof');

%evaluate if Clinical or Preclinical
if file(2) ==  51 && file(8) == [21041]
    %Clinical
    modality = 'Clinical';
    pSize = 51;
    channels = 72;
    text_header = 21041;
elseif file(2) ==  93 && file(8) == [21040]
    %Preclinical
    modality = 'Preclinical';
    pSize = 93;
    channels = 36;
    text_header = 21040;
else
    modality = 'Error';
end
disp(strcat('modality=',num2str(modality)));

%estimated number of events, considering 16bits words
estimated_n_events=floor((pSize*FileSize/2)/((8+channels)*pSize+4));
if nargin == 3
    estimated_n_events = num_events;
    file_pos = FileSize -2*(estimated_n_events*(8+channels)+4);
else
    file_pos = 0;
end
fseek(fid,file_pos,'bof');

%% Data Loading

count = 1;
Data = zeros(estimated_n_events,channels);
nodeID = zeros(estimated_n_events,1);
timestamp = zeros(estimated_n_events,1);
data_index = [];
for i =1:pSize
    data_index = [data_index (i-1)*(8+channels)+[9:8+channels]]; %event data peak matrix
end
heads_index = 5:8+channels:pSize*(8+channels); %nodeID
time_index1 = heads_index+1;
time_index2 = heads_index+2;
time_index3 = heads_index+3;

tic
while  ftell(fid)<FileSize 
    
    Header0=fread(fid,4,'uint16','b'); %data whole packet header
    if isequal(Header0,[0;pSize;21321;17495])  %0|pSize|SI|DW
        file_temp=fread(fid,(8+channels)*pSize,'uint16','b');
        nodeID(count:count+pSize-1)=file_temp(heads_index) / 8 / 256; %5 bits node ID 
        timestamp(count:count+pSize-1)=time_res*(file_temp(time_index1).*2^32+file_temp(time_index2).*2^16+file_temp(time_index3));
        vector=file_temp(data_index);
        Data(count:count+pSize-1,:)=reshape(vector,channels,pSize)';
        count = count+pSize;
    end
    
end

%delete lines with no real meaning
Data(nodeID==0,:)=[];
timestamp(nodeID==0)=[];
nodeID(nodeID==0)=[];

%ADC to ASIC reorder
[~,ASIC_order]=INSERT_reorder(modality);
Data(:,ASIC_order)=Data;
toc
fclose(fid);
end