function [output]=PlotSpectra(Frame,num_nodes,Node,bins,output,satur_level)

h1 = figure('units','normalized','outerposition',[0 0 1 1]);
h2 = figure('units','normalized','outerposition',[0 0 1 1]);

if num_nodes <= 1
    disp('Energy spectrum')
    ht = suptitle('Signal Spectrum');
else
    disp('Energy spectra')
    ht = suptitle('Signal Spectra');
end
set(ht,'FontSize',18,'FontWeight','bold')
cols = ceil(sqrt(num_nodes));
rows= ceil(num_nodes/cols);
for n = 1:num_nodes
    figure(h1)
    subplot(rows,cols,n)
    Frame_node = Frame(Node == n,:);
    %%%%%%%%%%%%%%%
    SIGNAL = sum(Frame_node,2);
    %%%%%%%%%%%%%%%
    [y_en, x_en] = hist(SIGNAL,bins);
    plot(x_en,y_en,'Linewidth',2)
    hold on
    SIGNAL=sum(Frame_node(sum(Frame_node>satur_level,2)>0,:),2); %sum(Frame_node>satur_level,2)>0 -> n° di ch che saturano>0 (basta un ch che satura) 
    [y_en_sat, x_en_sat] = hist(SIGNAL,bins);
    plot(x_en_sat,y_en_sat,'r','Linewidth',2)
    legend(['Node ',num2str(n)],'Saturations')
    xlabel('Signal sum [ADC channels]')
    ylabel('Counts [a.u.]')
    %save data into output variable
    output(n).spectrum.x = x_en;
    output(n).spectrum.y = y_en;
    
    figure(h2);
    subplot(rows,cols,n)
    SIGNAL=sum(Frame_node(sum(Frame_node>satur_level,2)==0,:),2); %eventi per cui nessun ch satura
    [y_en, x_en] = hist(SIGNAL,bins);
    plot(x_en,y_en,'Linewidth',2)
    legend(['Node ',num2str(n)])
    xlabel('Signal sum [ADC channels]')
    ylabel('Counts [a.u.]')
    %save data into output variable
    output(n).spectrum_filt.x = x_en;
    output(n).spectrum_filt.y = y_en;
end
end