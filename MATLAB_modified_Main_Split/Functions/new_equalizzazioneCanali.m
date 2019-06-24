function [posizioni,coefficienti]=new_equalizzazioneCanali(impulsazione,Frame)

if isrow(impulsazione)
    impulsazione = impulsazione';
end

numch = min(size(Frame));
numpeaks = numel(impulsazione) * numch;
lunghezza = max(max(Frame))+1;
allHist = zeros(numch*lunghezza, 1);
interval = 1:lunghezza;
% interval=1:2^12;
for i=1:numch
    allHist((i-1)*lunghezza + interval) = hist(Frame(:,i), interval);    
end

%filtro e tolgo la baseline (eventi di intensità < al 40% di tutti gli eventi)
% allHist = medfilt1(allHist,10);
% e1 = sort(allHist);
% e2 = cumsum(e1);
% e2 = e2/e2(end);
% allHist(allHist <  e1(dsearchn(e2,0.4))) = 0;

%filtraggio a media mobile
% windowSize = 8; %5
% b = (1/windowSize)*ones(1,windowSize);
% allHist = filter(b,1,allHist);

[value, idx, width, prominence] = findpeaks(allHist,'MinPeakDistance',100);

[~, width_idx] = sort(width,'descend');
[~, prominence_idx] = sort(prominence,'descend');

% if numpeaks > numel(width_idx)
%     C = intersect(width_idx(1:numel(width_idx)), prominence_idx(1:numel(width_idx)));
% else
    C = intersect(width_idx(1:numpeaks), prominence_idx(1:numpeaks));
% end

%figure, hist(ceil(idx(C)/lunghezza), numch)
rrrr = hist(ceil(idx(C)/lunghezza), numch);

figure('units','normalized','outerposition',[0.05 0.35 0.9 0.5]), 
subplot(1,5,1)
hist(ceil(idx(C)/lunghezza), numch)
title('n° picchi vs channels');
xlabel('channel');
ylabel('n° picchi identificati');
axis tight
subplot(1,5,2:5)
hold on
plot(allHist)
scatter(idx(C), value(C), 'filled')
peakpoints = idx(C);
title('valori in [ADCchannel] per ogni channel, in fila');
for i=1:numch
    if (rrrr(i) < numpeaks/numch)
        title(['in channel ', num2str(i), ' there are ', num2str(rrrr(i)), ' peaks /', num2str(numpeaks/numch)])
        
        xlim([ (i-1)*lunghezza+1 , i*lunghezza ])
        ylim([ 0 , 1.25*max(allHist) ])
        xpoints = ((i-1)*lunghezza+1 : i*lunghezza)';
                
        flag_continue = 0;
        while (flag_continue == 0)
            choise = menu('execute:', 'continue', 'select missing peak', 'exit');
            if (choise == 1)
                flag_continue = 1;
            end
            if (choise == 2)
                [xpoint, ypoint] = ginput(1);
                [fitting,~] = FcnGaussFitV2(xpoints, allHist(xpoints), round(xpoint), round(ypoint), 0.7, 0.7 );
                scatter(round(fitting.b1), fitting.a1, 'filled', 'red')
                legend off
                peakpoints = [peakpoints; round(fitting.b1)];
            end
            if (choise == 3)
                return
            end
        end   
    end
end

figure, hold on
plot(allHist)
scatter(peakpoints, allHist(peakpoints), 'filled')

%% rimozione picchi
% es. 8 picchi di impulsazione più 1 doppio picco da eliminare

for i=1:numch
    if (rrrr(i) > numpeaks/numch)
        rrrr = hist(ceil(peakpoints/lunghezza), numch);
        figure, hist(ceil(peakpoints/lunghezza), numch);

%         xpoints = ((i-1)*lunghezza+1 : i*lunghezza)';
%         figure,
%         plot(xpoints, allHist(xpoints))

%         gigio = peakpoints(logical(((peakpoints > (i-1)*lunghezza+1) .* (peakpoints < i*lunghezza))));
%         
%         prompt = {'Enter matrix size:','Enter colormap name:'};
%         dlg_title = 'Input';
%         num_lines = 1;
%         defaultans = {'20'};
%         answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        
    end
end

%% posizioni
allHist_x = repmat(interval, 1, numch)';
allHist_c = ceil((1:numch*lunghezza)/lunghezza)';

posizioni=[allHist_x(peakpoints),allHist_c(peakpoints)];

%% coefficienti (m,q)
k = 0; M = 0; Q = 0; problemi = 0; probCh = [];
for i=1:numch
    points = posizioni(posizioni(:,2) == i, 1);
    if (numel(points) == numpeaks / numch)
        fitting = fit(impulsazione, points, 'poly1');
        coefficienti.(['canale_',num2str(i)]).m = fitting.p1;
        coefficienti.(['canale_',num2str(i)]).q = fitting.p2;
        M = M + fitting.p1;
        Q = Q + fitting.p2;
        k = k + 1;
    else
        problemi = 1;
        probCh = [probCh; i];
    end
end
M = M/k;
Q = Q/k;

if (problemi)
    disp('problemi');
    for i=probCh
        value = impulsazione*M + Q;
        problematicPosition = posizioni(posizioni(:,2) == i, 1);


        % caso A: ci sono più valori di impulsazione che picchi
        [~,d] = dsearchn(value, problematicPosition);
        [s,~] = dsearchn(value, problematicPosition - min(d));
        fitting = fit(impulsazione(s), problematicPosition, 'poly1');
        coefficienti.(['canale_',num2str(i)]).m = fitting.p1;
        coefficienti.(['canale_',num2str(i)]).q = fitting.p2;
        figure, hold on
        scatter(impulsazione*fitting.p1+fitting.p2, ones(numpeaks / numch, 1));
        scatter(problematicPosition, ones(numel(problematicPosition), 1));
        
        % caso B: ci sono più picchi che valori di impulsazione
        % caso possibile se riconosce il rumore come picco
        % tuttavia se trovasse 8+1 picchi, avrei già eliminato il picco spurio sopra
        % credo pertanto sia irrealizzabile e dunque non implementato
        % ad ogni modo, dovrebbe essere lo 'speculare' del caso A

        
    end
end








end

