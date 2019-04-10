%% bulma RAW data path:
mainDir='C:\Users\INSERT\Documents\insert\matlab\quick\Dataset\20171002\';
%FileName = 'gamma_asics_collimator_node10_hv34v8_gain10th25_47mbq_10h28_planar_final_02';
%FileName = 'gamma_asics_collimator_nodeall_hv34v8_gain10th25_50mbq_10h00_uniformPhantom_11';
%FileName = 'gamma_asics_collimator_nodeall_hv34v8_gain10th25_14mbq9_12h42_2dHoffmanPhantom_04';
%FileName = 'gamma_asics_collimator_nodeall_hv34v8_gain10th25_37mbq4_15h45_JaszczakPhantom_09';
FileName = 'gamma_asics_collimator_nodeall_hv34v8_gain10th25_2halfLineSources_4mbq7_7mbq2_15h56_transaxCalib_final_04';
%FileName = 'bulmaraw_H01_F';
%FileName = 'gamma_asics_node01_hv34v8_gain10th25_4mbq_x_01';
FileExtension= '.data';
bulmarawFilePath=[mainDir , FileName , FileExtension ];

%% Target Directory for the converted digital raw
digitalrawTargetDir = ['C:\Users\INSERT\Documents\Clinical_insert_calibgui\_AUTOPlay\' , FileName ];
mkdir(digitalrawTargetDir)

%% Number of nodes
% always:
NumofNodes=10;

%% CONVERTing
MexPlayDataConverter(bulmarawFilePath,digitalrawTargetDir,NumofNodes);

%% Image saving File Path
saveimageFilePath = [digitalrawTargetDir , '\images.mat'];

%% Directory of the corr datas
CorrDir = 'C:\Users\INSERT\Documents\Clinical_insert_calibgui\_AUTOPlay\Corr' ;

%% PLAYing
%image = Fnc_PlayAcquisition(CorrDir,digitalrawTargetDir,saveimageFilePath);

image = Fnc_PlayAcquisition_n(CorrDir,digitalrawTargetDir,saveimageFilePath,29);
