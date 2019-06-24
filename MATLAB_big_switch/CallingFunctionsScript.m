%% All Possible Operations
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
% 44       Get Calibration Lines

%% Useful Folders Pathways
addpath('C:\Users\Lab FIORINI\Desktop\new insert\Functions')
addpath('C:\Users\Lab FIORINI\Desktop\new insert\Geometries')

%% Default Pathmanes
if ~exist('DefaultPathnames','var')
    DefaultPathnames.default = pwd;
    DefaultPathnames.DataPathname = DefaultPathnames.default;
    DefaultPathnames.CristalPathname = DefaultPathnames.default;
    DefaultPathnames.ChannelPathname = DefaultPathnames.default;
    DefaultPathnames.EnergyPathname = DefaultPathnames.default;
    DefaultPathnames.Ba133Pathname = DefaultPathnames.default;
    DefaultPathnames.OutputPathname = [pwd, '\Outputs\'];
end

%% Flags

% NUMERO DI EVENTI comment to see all the events
clear num_events
num_events = 5*10^6;

% TILE
flag.Tile = 0;

% BASIC
flag.histograms = 1;
flag.spectra = 1;
flag.images = 1;
flag.baseline = 600;
flag.multipleEW = 0;

% CALIBRAZIONE CRISTALLO
flag.equalization_cristallo = 0;

% CALIBRAZIONE CANALI
flag.equalization = 1;
flag.equalization_show = 1;
flag.equalization_use = 1;

% SPETTRO LOCALE
flag.local_spectrum = 1;
flag.sector = 3;

% % CALIBRAZIONE IN ENERGIA
% flag_EgCalibration = 0;
% flag_EgCalibration_use = 1;

% GAUSSIAN FIT E RISOLUZIONE ENERGETICA IN CANALI ADC
flag.RS = 1;

% CALIBRAZIONE ENERGETICA
% carica output del Ba133 precedentemente elaborato per calcolare la calibrazione energetica
flag.CE = 1;

% FRAME FILTERING
flag.step = 100; flag.frameFilteringRE = 0;

% ALTRO

%% Save in ExecutionDiary.txt the stdout
diary off
delete ExecutionDiary.txt
diary('ExecutionDiary.txt')
diary on
disp('---------------STARTING POINT---------------')
disp(datetime('now'))

%% Initialize Variables
output = [];
ws = [];
if exist('num_events','var')
    ws.num_events = num_events;
else
    ws.num_events = NaN;
end
ws.flag = flag;
ws.DefaultPathnames = DefaultPathnames;
[ output,ws ] = do_fcn_type( 'clear', output,ws );


%% other variables
load('spots_ROI_Grid_2601.mat')
ws.spots_ROI = spots_ROI;