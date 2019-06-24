function [spectrum,spectrum_filt]=PlotSpectrum(Frame,NODE,bins,satur_level)
disp('Energy spectrum')

figure('units','normalized','outerposition',[0 0 1 1]);
%%%%%%%%%%%%%%%
SIGNAL = sum(Frame,2);
%%%%%%%%%%%%%%%
[y_en, x_en] = hist(SIGNAL,bins);
plot(x_en,y_en,'Linewidth',2)
hold on
SIGNAL=sum(Frame(sum(Frame>satur_level,2)>0,:),2); %sum(Frame>satur_level,2)>0 -> n° di ch che saturano>0 (basta un ch che satura) 
[y_en_sat, x_en_sat] = hist(SIGNAL,bins);
plot(x_en_sat,y_en_sat,'r','Linewidth',2)
legend(['Node ',num2str(NODE)],'Saturations')
xlabel('Signal sum [ADC channels]')
ylabel('Counts [a.u.]')
title('Signal Spectrum')
set(gca,'FontSize',18,'FontWeight','bold')
%save data into output variable
spectrum.x = x_en;
spectrum.y = y_en;
    
figure('units','normalized','outerposition',[0 0 1 1]);
%%%%%%%%%%%%%%%
SIGNAL=sum(Frame(sum(Frame>satur_level,2)==0,:),2); %eventi per cui nessun ch satura
%%%%%%%%%%%%%%%
[y_en, x_en] = hist(SIGNAL,bins);
plot(x_en,y_en,'Linewidth',2)
legend(['Node ',num2str(NODE)])
xlabel('Signal sum [ADC channels]')
ylabel('Counts [a.u.]')
title('Signal Spectrum without Saturation')
set(gca,'FontSize',18,'FontWeight','bold')
%save data into output variable
spectrum_filt.x = x_en;
spectrum_filt.y = y_en;
end