function [ fitting ] = monoG_v2( bins, spettro, xpoint, ypoint )

if isrow(spettro)
    spettro = spettro';
end
if isrow(bins)
    bins = bins';
end

% figure('units','normalized','outerposition',[0 0 1 1],'Name','ginput')
% plot(bins, spettro)
% pause(wait_time)
% [q,w] = ginput(1);              % q = posizione
% close ginput

Soglia_fit_dx = 0.8;
Soglia_fit_sx = 0.6;
%y_flt = medfilt1(spettro,7);
y_flt = spettro;
x_max = round(xpoint);
y_max = round(ypoint);
[fitting,~] = FcnGaussFitV2(bins,y_flt,x_max,y_max,Soglia_fit_dx,Soglia_fit_sx );

% figure('units','normalized','outerposition',[0 0 1 1],'Name','ginput')
% plot(fitting, bins, spettro)

cla
plot(fitting, bins, spettro)

end

