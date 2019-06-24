function [Data,Head_index,Time_stamp,modality] = openDataFile(fid,num_events)
clc
disp('Loading .data File')
%initializations and parameters
time_res = 25e-9;
%LUT = [0;0;0;0;0;0;0;1;11;0;0;0;0;0;0;2;12;0;0;0;0;0;0;3;13;0;0;0;0;0;0;4;14;0;0;0;0;0;0;5;15;0;0;0;0;0;0;6;16;0;0;0;0;0;0;7;17;0;0;0;0;0;0;8;18;0;0;0;0;0;0;9;19;0;0;0;0;0;0;10;20];


%% File opening
%fid = fopen([pwd,File_Name],'r','b'); % windows
%check file dimension
fseek(fid,0,'eof');
FileSize=ftell(fid);
fseek(fid,0,'bof');

file = fread(fid,8,'uint16=>uint16','b');
fseek(fid,0,'bof');

% evaluate if Clinical or Preclinical
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
% Estimated number of events, considering 16bits words
estimated_n_events =floor( pSize*FileSize/2/((8+channels)*pSize+4) );
if nargin == 2
    estimated_n_events = num_events;
    file_pos = FileSize -2*(estimated_n_events*(8+channels)+4);
else
    file_pos = 0;
end
fseek(fid,file_pos,'bof');
%% Data Loading
Count = 1;
Data = zeros(estimated_n_events,channels);
Head_index = zeros(estimated_n_events,1);
Time_stamp = zeros(estimated_n_events,1);
data_index = [];
for i =1:pSize
    data_index = [data_index (i-1)*(8+channels)+[9:8+channels]];
end
heads_index = 5:8+channels:pSize*(8+channels);
time_index1 = heads_index+1;
time_index2 = heads_index+2;
time_index3 = heads_index+3;
tic
while  ftell(fid)<FileSize 
    
    Header0=fread(fid,4,'uint16','b');
    if isequal(Header0,[0;pSize;21321;17495])
        file_temp=fread(fid,(8+channels)*pSize,'uint16','b');
        Head_index(Count:Count+pSize-1)=file_temp(heads_index) / 8 / 256;
        %Head_index(Count:Count+pSize-1)=LUT(file_temp(heads_index)/256);
        Time_stamp(Count:Count+pSize-1)=time_res*(file_temp(time_index1).*2^32+file_temp(time_index2).*2^16+file_temp(time_index3));
        Vector=file_temp(data_index);
        Data(Count:Count+pSize-1,:)=reshape(Vector,channels,pSize)';
        Count = Count+pSize;
    end
    
end
%delete lines with no real meaning
Data(Head_index==0,:)=[];
Time_stamp(Head_index==0)=[];
Head_index(Head_index==0)=[];

%ADC to ASIC reorder
 [~,ASIC_order]=INSERT_reorder(modality);
 Data(:,ASIC_order)=Data;
toc
fclose(fid);
end