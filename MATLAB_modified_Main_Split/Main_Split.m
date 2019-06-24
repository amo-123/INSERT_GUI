%% Main to Split Nodes
%Code to load .data files (INSERT files) and split the acquisitions from different nodes into FRAME_NODES variable

close all
clear all
clc

if ~exist('current_path', 'var')
    current_path = pwd;
end
current_path=strcat(current_path,'\');
addpath('Functions');
addpath('Geometries');
addpath(strcat('Functions\.data Functions'));

%% Initializations

%determine if the user wants to open only the last num_events
num_events=5000000; %number of events to display (comment this line to see all the events)

%% Load

FilterSpec = '*.data';
[filename,filepath] = uigetfile([pwd,'\',FilterSpec], 'Select .data file', 'MultiSelect', 'off');

%load function
if exist('num_events','var')
   [Frame,Node,Time_stamp,modality]=openDataFile(filename,filepath,num_events);
   disp(strcat('last',{' '},num2str(num_events),{' '},'events loaded'))
else
   [Frame,Node,Time_stamp,modality]=openDataFile(filename,filepath);
   disp('all events loaded')
end

%load proper reorder array and geometrical reorder
[Phys_order]=INSERT_reorder(modality); 
Frame(:,Phys_order)=Frame;
Node=round(Node); 
num_nodes=max(Node); %number of nodes

%divide datasets of different nodes
FRAME_NODE{num_nodes,1}=0;
for n=1:num_nodes  
    FRAME_NODE{n,1}=Frame(Node==n,:);
end

disp(strcat('Modality:',32,modality,32,'-> Number of Nodes:',32,num2str(num_nodes)))

%% Check number of events x node

for n=1:num_nodes
    Events_counts(n)=length(FRAME_NODE{n,1});
end
figure
plot(1:1:num_nodes,Events_counts,'-o','linewidth',2)
xlabel('#node')
ylabel('#events')
set(gca,'fontsize',15,'fontweight','bold')
grid on

Tot_events=sum(Events_counts);
