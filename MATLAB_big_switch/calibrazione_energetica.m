function [ cobalto, bario ] = calibrazione_energetica( output, pathname_Ba133, filename_Ba133 )

%% ex pamela
%% colcolo della risoluzione energetica del Co57 calibrato col Ba133
%% entrambi devono essere stati elaborati dal main

%% variabili globali
fotopicco_Co57 = 122; % keV
fotopicco_Ba133 = 81; % keV

%% carico Co57
cobalto = output;

%% carico Ba133
bario = load([pathname_Ba133,filename_Ba133]);
bario = bario.(cell2mat(fieldnames(bario)));

%% calcolo retta calibrazione ADC -> EG
m = (fotopicco_Co57 - fotopicco_Ba133) / (cobalto.peak_fitting{1,1}.b1 - bario.peak_fitting{1,1}.b1);
q = fotopicco_Co57 - m*cobalto.peak_fitting{1,1}.b1;

%% calibrazione
bins = cobalto.spectrum.x;
energy_axis = 0:0.05:204.8;

x_EG = m*bins + q;

spectrum_coba = interp1(x_EG, cobalto.spectrum.y, energy_axis, 'linear');
spectrum_baro = interp1(x_EG, bario.spectrum.y, energy_axis, 'linear');

figure, axis('tight'), title('Spettri calibrati in energia di Co57 e Ba133'), hold on
plot(energy_axis, spectrum_coba/max(spectrum_coba))
plot(energy_axis, spectrum_baro/max(spectrum_baro))

%% calcolo della risoluzione energetica
bario.peak_fitting_calibrato{1} = monoG(energy_axis, spectrum_baro, 0);
bario.peak_fitting_calibrato{2} = 2.3548*bario.peak_fitting_calibrato{1}.c1/sqrt(2)/bario.peak_fitting_calibrato{1}.b1;

cobalto.peak_fitting_calibrato{1} = monoG(energy_axis, spectrum_coba, 0);
cobalto.peak_fitting_calibrato{2} = 2.3548*cobalto.peak_fitting_calibrato{1}.c1/sqrt(2)/cobalto.peak_fitting_calibrato{1}.b1;

title(['la risoluzione energetica del fotopicco Co57 è ', num2str(100*cobalto.peak_fitting_calibrato{2}), ' %'])


end

