function [output]=PlotSectra(Frame,num_nodes,Node,bins,FName)
h1 = figure('units','normalized','outerposition',[0 0 1 1],'Name',FName,'NumberTitle','off');
h2 = figure('units','normalized','outerposition',[0 0 1 1],'Name',FName,'NumberTitle','off');
% satur_level = 2^12-20;
satur_level = 4094;

clc
if num_nodes <= 1
    disp('Energy spectrum')
    ht = suptitle('Signal Spectrum');
else
    disp('Energy spectra')
    ht = suptitle('Signal Spectra');
end
set(ht,'FontSize',18,'FontWeight','bold');
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
    SIGNAL =sum(Frame_node(sum(Frame_node>satur_level,2)>0,:),2);
    [y_en_sat, x_en_sat] = hist(SIGNAL,bins);
    plot(x_en_sat,y_en_sat,'r','Linewidth',2)
    legend(['Node ',num2str(n)],'Saturations')
    xlabel('Signal sum [ADC channels]')
    ylabel('Counts [a.u.]')
    %save data into output variable
    output.(['node_',num2str(n)]).x = x_en;
    output.(['node_',num2str(n)]).y = y_en;
    
    figure(h2)
    subplot(rows,cols,n)
    SIGNAL =sum(Frame_node(sum(Frame_node>satur_level,2)==0,:),2);
    [y_en, x_en] = hist(SIGNAL,bins);
    plot(x_en,y_en,'Linewidth',2)
    legend(['Node ',num2str(n)])
    xlabel('Signal sum [ADC channels]')
    ylabel('Counts [a.u.]')
end


end