function [ num_fnc_type ] = conversione_fcnType2num( string_fcn_type )

corrispondence_num(1) = 21;
corrispondence_str{1} = 'clear';

corrispondence_num(2) = 20;
corrispondence_str{2} = 'select and load .data file';

corrispondence_num(3) = 1;
corrispondence_str{3} = 'load .data';

corrispondence_num(4) = 2;
corrispondence_str{4} = 'num nodi e riordino .data';

corrispondence_num(5) = 3;
corrispondence_str{5} = 'divido nei vari nodi';

corrispondence_num(6) = 4;
corrispondence_str{6} = 'Equalizatione del cristallo';

corrispondence_num(7) = 5;
corrispondence_str{7} = 'Equalizatione dei canali';

corrispondence_num(8) = 6;
corrispondence_str{8} = 'Histograms';

corrispondence_num(9) = 7;
corrispondence_str{9} = 'Signal Spectra';

corrispondence_num(10) = 8;
corrispondence_str{10} = 'Gaussian Fitting e Risoluzione Energetica';

corrispondence_num(11) = 9;
corrispondence_str{11} = 'Images (Modified Centroid Method Reconstruction)';

corrispondence_num(12) = 10;
corrispondence_str{12} = 'Spettri Locali';

corrispondence_num(13) = 11;
corrispondence_str{13} = 'CALIBRAZIONE ENERGETICA';

corrispondence_num(14) = 12;
corrispondence_str{14} = 'SAVE FRAME AS .mat FILE';

corrispondence_num(15) = 13;
corrispondence_str{15} = 'Save Output';

corrispondence_num(16) = 36;
corrispondence_str{16} = 'Local Spectra by Spots';

corrispondence_num(17) = 44;
corrispondence_str{17} = 'Get Calibration Lines';

corrispondence_num(18) = 51;
corrispondence_str{18} = 'frame filtering';

corrispondence_num(19) = 19;
corrispondence_str{19} = 'Spettri Locali 2';

corrispondence_num(20) = 29;
corrispondence_str{20} = 'sottrazione offset e calibrazione monopicco';

%% find num_fnc_type
num_fnc_type = corrispondence_num( strcmp(corrispondence_str, char(string_fcn_type)) );

if isempty(num_fnc_type)
    num_fnc_type = NaN;
end

end

% exemple calling function
% 21       clear
% 20       start: select .data file
% 1        load .data
% 2        num nodi e riordino .data
% 3        divido nei vari nodi
% 4        Equalizatione del cristallo
% 5        Equalizatione dei canali
% 6        Histograms
% 7        Signal Spectra
% 8        Gaussian Fitting e Risoluzione Energetica
% 9        Images (Modified Centroid Method Reconstruction)
% 10       Spettri Locali
% 11       CALIBRAZIONE ENERGETICA
% 12       SAVE FRAME AS .mat FILE
% 13       Save Output