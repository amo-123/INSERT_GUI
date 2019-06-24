%% unisci file per calibrazione canale
% close all
% clearvars
% clc
FilterSpec = '*data';
[DataFilename, pathname] = uigetfile(FilterSpec, 'select .data acquisition file', 'MultiSelect', 'on');

F = [];
if (iscell(DataFilename))
    for i=1:size(DataFilename,2)
        filename = DataFilename{1,i};
        fid = fopen([pathname,filename],'r','b'); % per win
        [Frame,Node,Time_stamp,modality]=openDataFile(fid);
        F = [F; Frame];
    end
else
    fid = fopen([pathname,DataFilename],'r','b'); % per win
    [Frame,Node,Time_stamp,modality]=openDataFile(fid);
    F = [F; Frame];
end

% clearvars -except F
% 
% [ ~, m_corr, q_corr ] = equalizz_ch( F, 1:72, 1 );
% 
% %save('20170920_calibImp_g12_th30_C03_7peaks.mat','F')
% 
% 
% nodo = 1;
% 
% b = b(1:72);
% m = m(1:72);
% 
% save(['q_moduloSingolo_',num2str(nodo),'.mat'], 'b')
% save(['m_moduloSingolo_',num2str(nodo),'.mat'], 'm')
% 
% clearvars
% close all
% 
% [ ~, m_corr, q_corr ] = equalizz_ch( F, 1:36, 1 );
