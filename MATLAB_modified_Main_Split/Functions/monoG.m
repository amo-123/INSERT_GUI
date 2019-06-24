function [fitting]=monoG(bins_spettro,counts_spettro,img_title,soglia_fit_dx,soglia_fit_sx)

interval=FcnUsefullData(counts_spettro,0);
counts_spettro = counts_spettro(interval);
bins_spettro = bins_spettro(interval);

if isrow(counts_spettro)
    counts_spettro = counts_spettro';
end
if isrow(bins_spettro)
    bins_spettro = bins_spettro';
end

figure('units','normalized','outerposition',[0 0 1 1],'Name','ginput')
plot(bins_spettro,counts_spettro,'linewidth',2)
title(img_title)
set(gca,'fontweight','bold','fontsize',18)
chioce=0;
while chioce==0
    chioce=menu('Zoom the peak and press OK before selecting the peak','OK');
end
[x_peak,~]=ginput(1); 
x_peak=round(x_peak);
idx_peak=find(bins_spettro>=x_peak,1);
if abs(bins_spettro(idx_peak-1)-x_peak)<abs(bins_spettro(idx_peak)-x_peak)
    idx_peak=idx_peak-1;
end
% x_peak=bins_spettro(idx_peak);
y_peak=counts_spettro(idx_peak);

[fitting,~]=FcnGaussFitV2(bins_spettro,counts_spettro,x_peak,y_peak,soglia_fit_dx,soglia_fit_sx);
%peak update 
x_peak=round(fitting.b1);
diff=abs(bins_spettro-x_peak);
[~,idx_diff]=sort(diff,'ascend');
idx_peak=idx_diff(1);
y_peak=counts_spettro(idx_peak);

hold on
plot(x_peak,y_peak,'xr','markersize',10,'linewidth',2)
plot(bins_spettro,fitting(bins_spettro),'--k','linewidth',2)
hold off

end

