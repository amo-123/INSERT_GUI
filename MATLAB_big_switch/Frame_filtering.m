function [ Frame_others ] = Frame_filtering(Frame, varargin)

small = 0;
re = 0;

thr_low = 10;
thr_high = 4000;

switch nargin
    case 1
        step = 20;
    case 2    
        step = varargin{1};
    case 3
        step = varargin{1};
        re = varargin{2};
    case 4
        step = varargin{1};
        re = varargin{2};
        if strcmp(varargin{end}, 'small')
            small = 1;
        end
end


% finestra di energie
[all, x_all] = hist(sum(Frame,2),1000);

if small
    xpoint = [x_all(1), x_all(end)];
else
    figure,
    plot(x_all, all)
    while( menu('Choise the Spectra figure, zoom and press the "ok" button below, then select the photopeak WINDOWS', 'ok') ~= 1)
    end
    [xpoint, ~] = ginput(2);
end

Frame = Frame( and(sum(Frame,2) > xpoint(1), sum(Frame,2) < xpoint(2)), : );

% scelgo numero di bins
[all, x_all] = hist(sum(Frame,2),1000);
all = cumsum(all);
bins = xpoint(1) : step : xpoint(2);

% trovo porzioni di Frame corrispondenti a zeri e saturazioni
[under10, x_under10] = hist(sum(Frame( min(Frame, [], 2) <= thr_low, : ),2),bins);
[over4095, x_over4095] = hist(sum(Frame( max(Frame, [], 2) >= thr_high, : ),2),bins);
Frame_others = Frame( and(min(Frame, [], 2) > thr_low, max(Frame, [], 2) < thr_high), : );
[others, x_others] = hist(sum(Frame_others,2),bins);
[all, x_all] = hist(sum(Frame,2),bins);
tot = sum(all)/100;

disp([['zeros = ',num2str(round(sum(under10)/tot,2)),' %'], 32, ['saurations = ', num2str(round(sum(over4095)/tot,2)),' %'], 32, ['num events = ', num2str(size(Frame_others,1))]])

if (re==1)
    %% risoluzione energetica

    % others
    plot(x_others, others)
    while( menu('Choise the Spectra figure, zoom and press the "ok" button below, then select the photopeak apex', 'ok') ~= 1)
    end
    [xpoint, ypoint] = ginput(1);
    fitting = monoG_v2(x_others, others, xpoint, ypoint);
    EnergyRisolution_others = 2.3548*fitting.c1/sqrt(2)/fitting.b1;
    EG_others = ['Energy resolution = ', num2str(round(EnergyRisolution_others*100,2)), ' %'];

    % all
    plot(x_all, all)
    while( menu('Choise the Spectra figure, zoom and press the "ok" button below, then select the photopeak apex', 'ok') ~= 1)
    end
    [xpoint, ypoint] = ginput(1);
    fitting = monoG_v2(x_others, others, xpoint, ypoint);
    EnergyRisolution_all = 2.3548*fitting.c1/sqrt(2)/fitting.b1;
    EG_all = ['Energy resolution = ', num2str(round(EnergyRisolution_all*100,2)), ' %'];

    % confronto di risoluzioni energetiche
    cla, hold on
    plot(x_all, all, 'black')
    plot(x_others, others, 'red', 'linewidth', 2)
    plot(x_under10, under10)
    plot(x_over4095, over4095)

    legend(['complessivo',10,EG_all], ['others = ',num2str(round(sum(others)/tot,2)),' %',10,EG_others], ['zeros = ',num2str(round(sum(under10)/tot,2)),' %'], ['saurations = ', num2str(round(sum(over4095)/tot,2)),' %'])
end

end

