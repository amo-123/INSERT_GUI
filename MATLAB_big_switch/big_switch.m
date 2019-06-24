if ~exist('selected_fnc','var')
    CallingFunctionsScript; %---------------------------------------------------> External Script Showing All Executable Operations
    %% Desired Operations
    
%     ws.flag.small = 'small'; 
% [ output,ws ] = do_fcn_type( 'frame filtering', output,ws );
%     [ output,ws ] = do_fcn_type( 'Histograms', output,ws );

    [ output,ws ] = do_fcn_type( 'Signal Spectra', output,ws );

    [ output,ws ] = do_fcn_type( 'Images (Modified Centroid Method Reconstruction)', output,ws );

%     [ output,ws ] = do_fcn_type( 'Spettri Locali 2', output,ws );
   
%     [ output,ws ] = do_fcn_type( 'sottrazione offset e calibrazione monopicco', output,ws );

%     [ output,ws ] = do_fcn_type( 'Gaussian Fitting e Risoluzione Energetica', output,ws );
%     [ output,ws ] = do_fcn_type( 'Local Spectra by Spots', output,ws );

%     [ output,ws ] = do_fcn_type( 'Get Calibration Lines', output,ws );
%     [ output,ws ] = do_fcn_type( 'select and load .data file', output,ws );
%     [ output,ws ] = do_fcn_type( 'divido nei vari nodi', output,ws );
    
    %% End of the code
    selected_fnc = 0;
end

for doing_function = selected_fnc
    if exist('output','var')
        ws.idx_fnc = numel(output)+1;
    else
        ws.idx_fnc = 1;
    end
    
    switch doing_function
        case 0
            %% End of the code
            DefaultPathnames = ws.DefaultPathnames;
            clear selected_fnc doing_functions
            clc, disp(['That',39,'s all folks! THE END!']), diary off
        case 22 
            %% flags
            
        case 51
            %% frame filtering
            % rimozione eventi contenenti zeri (valori <= 10) e/o saturazioni (valori >= 4095)
            output{ws.idx_fnc}.fcn_type = 'frame filtering';
            [node_check, output,ws] = get_field('node_check', 'select and load .data file' ,output,ws);
            [ idxFrame ] = idx_done_fcn( 'select and load .data file', output );
            
            for n = node_check
                [Frame_node, output,ws] = get_field( ['Frame_node.node_',num2str(n)], 'divido nei vari nodi', output,ws,'ref' );
                [ idx ] = idx_done_fcn( 'divido nei vari nodi', output );
                disp(n)
                output{idx}.Frame_node.(['node_',num2str(n)]) = Frame_filtering(eval(Frame_node), ws.flag.step, ws.flag.frameFilteringRE, ws.flag.small);
            end
            
            if numel(node_check) == 1
                output{idxFrame}.Frame = output{idx}.Frame_node.(['node_',num2str(n)]);
                ws.Node = node_check * ones(size(output{idxFrame}.Frame,1),1);
            else
                output{idxFrame}.Frame = [];
                ws.Node = [];
                for n = node_check
                    output{idxFrame}.Frame = [output{idxFrame}.Frame; output{idx}.Frame_node.(['node_',num2str(n)])];
                    ws.Node = [ws.Node; n * ones(size(output{idx}.Frame_node.(['node_',num2str(n)]),1),1)];
                end
            end
        case 29
            %% sottrazione offset e calibrazione monopicco
            output{ws.idx_fnc}.fcn_type = 'sottrazione offset e calibrazione monopicco';
            [node_check, output,ws] = get_field('node_check', 'select and load .data file' ,output,ws);
            for n = node_check
                [spectrum.x, output,ws] = get_field( ['spectrum.node_',num2str(n),'.x'], 'Signal Spectra', output,ws );
                [Frame_node, output,ws] = get_field( ['Frame_node.node_',num2str(n)], 'divido nei vari nodi', output,ws );
                bins = length(spectrum.x);
                [ idx ] = idx_done_fcn( 'Signal Spectra', output );
                % seleziono file .mat contenente gli offset
                [FileName,PathName,~] = uigetfile('.mat',['select offset file for node ',num2str(n)]);
                offset = load([PathName,FileName]);
                offset = offset.(cell2mat(fieldnames(offset)));
                offset = offset(1:ws.ChNum);
                if iscolumn(offset)
                    offset = offset';
                end
                % sottraggo gli offset
                Frame_node = Frame_node - repmat(offset,size(Frame_node,1),1);
                % fitting monogaussiano del fotopicco
                [y,x] = hist(sum(Frame_node,2), bins);
                figCal = figure();
                plot(x,y)
                % inserisci parametri monoG ed energia fotopicco
                while( menu('Choise the Spectra figure, zoom and press the "ok" button below, then select the photopeak apex', 'ok') ~= 1)
                end
                prompt = {'Enter monoG Sx thr:','Enter monoG Dx thr:','Enter photopeak energy:'};
                dlg_title = 'Input';
                num_lines = 1;
                defaultans = {'0.3','0.8','122'};
                answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                [q,w] = ginput(1);
                %close gcf
                [fitting,~] = FcnGaussFitV2(x', y', round(q), round(w), str2double(answer{2}), str2double(answer{1}) );
                subplot(2,1,1), plot(fitting, x,y)
                % slope della retta di calibrazione energetica
                do_fitting = cfit(fittype('poly1'), str2double(answer{3}) / fitting.b1, 0);
                x = do_fitting(x);
                subplot(2,1,2), plot(x,y)
                % salvo in spectrum
                output{idx}.spectrum.(['node_',num2str(n)]).x = x;
            end

        case 21
            %% clear
            close all
            clearvars -except DefaultPathnames ws
            clc
            addpath('Functions')
            addpath('Geometries')
            ws.idx_fnc = 1;
            %% output
            output{ws.idx_fnc}.fcn_type = 'clear';
        case 20
            %% start: select .data file
            FilterSpec = '*.data';
            [DataFilename, ws.pathname] = uigetfile([ws.DefaultPathnames.DataPathname,'\',FilterSpec], 'select .data acquisition file', 'MultiSelect', 'on');
            ws.DefaultPathnames.DataPathname = ws.pathname;
            ws.filename = DataFilename;
            %% output
            output{ws.idx_fnc}.fcn_type = 'select and load .data file';
            %% load .data
            fid = fopen([ws.pathname,ws.filename],'r','b'); % per win
            ws.FName = ws.filename(1:end-length(FilterSpec)+1);
            disp(['ACQUISITION FILE at ',ws.pathname,ws.filename])
            if ~isnan(ws.num_events)
                [Frame,Node,Time_stamp,modality]=openDataFile(fid,ws.num_events);
            else
                [Frame,Node,Time_stamp,modality]=openDataFile(fid); 
            end                
            %% num nodi e riordino .data
            ws.Node = round(Node);
            ws.num_nodes = max(Node);
            node_check = find(ismember(1:max(Node), Node));
            ws.node_check = node_check;
            [Phys_order]=INSERT_reorder(modality);
            Frame(:,Phys_order) = Frame;
            %ws.Frame = Frame;
            ChNum = numel(Phys_order);
            ws.ChNum = ChNum;
            %% save general information
            output{ws.idx_fnc}.fcn_type = 'select and load .data file';
            output{ws.idx_fnc}.Frame = Frame;
            output{ws.idx_fnc}.pathname = ws.pathname;
            output{ws.idx_fnc}.filename = ws.filename;
            output{ws.idx_fnc}.FName = ws.FName;
            output{ws.idx_fnc}.modality = modality;
            ws.modality = modality;
            output{ws.idx_fnc}.node_check = node_check;
            output{ws.idx_fnc}.Phys_order = Phys_order;
            ws.Phys_order = Phys_order;
            output{ws.idx_fnc}.ChNum = ChNum;
        case 1
            
        case 101
            %% save case1
            output{ws.idx_fnc}.fcn_type = 'save case1';
            output{ws.idx_fnc}.Frame = Frame;
            output{ws.idx_fnc}.Node = Node;
        case 3
            %% divido nei vari nodi
            output{ws.idx_fnc}.fcn_type = 'divido nei vari nodi';
            FrameRef = get_field( 'Frame', 'select and load .data file', output,ws,'ref' );
            [ idx ] = idx_done_fcn( 'divido nei vari nodi', output );
            for n=1:ws.num_nodes
                output{idx}.Frame_node.(['node_', num2str(n)]) = eval([FrameRef, '(ws.Node == ', num2str(n), ',:);']);
                % output{idx}.Frame_node.(['node_', num2str(n)]) = Frame(ws.Node == n,:);
                % Frame_node = Frame(Node == n,:);
                % nome = [FName, '_', num2str(n)];
                % save(nome, 'Frame_node')
            end
        case 4
            %% Equalizatione del cristallo
            if (flag.equalization_cristallo == 1)
                clc 
                disp('Channels Calibration')
                % carico o calcolo i coefficienti
                FilterSpec = '*.mat';
                [Filename, pathname] = uigetfile([DefaultPathnames.CristalPathname,'\',FilterSpec], ['select .mat CRISTAL equalization file for node from 1 to ',num2str(num_nodes)], 'MultiSelect', 'off');
                DefaultPathnames.CristalPathname = pathname;
                ChEq = load([pathname,Filename]);
                ChEq = ChEq.(cell2mat(fieldnames(ChEq)));
                for n=node_check
                    output(n).CristalCorrections = ChEq.(['node_',num2str(n)]);
                end
                % corrego il Frame
                for n=node_check
                    m_corr = output(n).CristalCorrections(:,1);
                    q_corr = output(n).CristalCorrections(:,2);

                    Frame_corr = Frame(Node == n,:);
                    Frame_corr(Frame_corr == 0) = NaN;
                    Frame_corr = Frame_corr + (rand(size(Frame_corr)) - 0.5);
                    for i=1:ChNum
                        Frame_corr(:,i)=Frame_corr(:,i)*m_corr(i) + q_corr(i);
                    end
                    Frame_corr = round(Frame_corr);
                    Frame_corr(Frame_corr < 0) = 0;
                    Frame_corr(Frame_corr > 4095) = 4095;
                    Frame_corr(isnan(Frame_corr)) = 0;
                    % correggo Frame originale
                    Frame(Node == n,:) = Frame_corr;
                end     
            else
                for n = node_check
                   output(n).CristalCorrections = [];
                end
            end
        case 5
            %% Equalizatione dei canali
            if (flag.equalization == 1)
                clc 
                disp('Channels Calibration')
                % carico o calcolo i coefficienti
                if (flag.equalization_use == 1)
                    FilterSpec = '*.mat';
                    [Filename, pathname] = uigetfile([DefaultPathnames.ChannelPathname,'\',FilterSpec], ['select .mat channels equalization file for node from 1 to ',num2str(num_nodes)], 'MultiSelect', 'off');
                    DefaultPathnames.ChannelPathname = pathname;
                    ChEq = load([pathname,Filename]);
                    ChEq = ChEq.(cell2mat(fieldnames(ChEq)));
                    for n=node_check
                        output(n).ChannelsCorrections = ChEq.(['node_',num2str(n)]);
                    end
                else
                    for n=node_check
                        Frame_corr = Frame(Node == n,:);
                        % caricamento
                        FilterSpec = '*.mat';
                        [Filename, pathname] = uigetfile(FilterSpec, ['select .mat channels calibration (chs x events) file for node ',num2str(n)], 'MultiSelect', 'off');
                        calibration = load([pathname,Filename]);
                        calibration = calibration.(cell2mat(fieldnames(calibration)));
                        % riordino della calibrazione
                        calibration(:,Phys_order) = calibration;
                        % calibrazione
                        [ ~, m_corr, q_corr ] = equalizz_ch( calibration, Phys_order, flag.equalization_show );
            %             ChEq = [m_corr, q_corr];
            %             output(n).ChannelsCorrections = ChEq;
            %             save(['ChEq_',Filename(1:end-4)], 'ChEq');
            %             

                        ChEq.(['node_',num2str(n)]) = [m_corr, q_corr];
                        output(n).ChannelsCorrections = ChEq.(['node_',num2str(n)]);
                        if (n==num_nodes)
                            save(['ChEq_','node_1-',num2str(n),'_',FName,'.mat'], 'ChEq');
                        end
                    end
                end
                % corrego il Frame
                for n=node_check
                    m_corr = output(n).ChannelsCorrections(:,1);
                    q_corr = output(n).ChannelsCorrections(:,2);

                    Frame_corr = Frame(Node == n,:);
                    Frame_corr(Frame_corr == 0) = NaN;
                    Frame_corr = Frame_corr + (rand(size(Frame_corr)) - 0.5);
                    for i=1:ChNum
                        Frame_corr(:,i)=Frame_corr(:,i)*m_corr(i) + q_corr(i);
                    end
                    Frame_corr = round(Frame_corr);
                    Frame_corr(Frame_corr < 0) = 0;
                    Frame_corr(Frame_corr > 4095) = 4095;
                    Frame_corr(isnan(Frame_corr)) = 0;
                    % correggo Frame originale
                    Frame(Node == n,:) = Frame_corr;
                end     
            else
                for n = node_check
                   output(n).ChannelsCorrections = [];
                end
            end
        case 6
            %% Histograms
            output{ws.idx_fnc}.fcn_type = 'histogram';
            if ws.flag.histograms == 1
                clc 
                disp('Single channels histograms')
                
                
                [node_check, output,ws] = get_field('node_check', 'select and load .data file' ,output,ws);
                f = figure('units','normalized','outerposition',[0 0 1 1],'Name',ws.FName,'NumberTitle','off');
              
                tab_group = uitabgroup('Parent', f);
                
                for n = node_check
                    [Frame_node, output,ws] = get_field( ['Frame_node.node_',num2str(n)], 'divido nei vari nodi', output,ws );
                    Frame_node = Frame_node + rand(size(Frame_node)) - 0.5;
                    
                    try
                        bins = linspace(0,round(max(max(Frame_node))),2^12);
                        output{ws.idx_fnc}.histograms.(['node_',num2str(n)])=PlotHistograms(ws.modality,Frame_node,ws.Phys_order,n,ws.FName,bins,tab_group,f);
                        drawnow
                    catch 
                        output{ws.idx_fnc}.histograms.(['node_',num2str(n)])=[];
                    end
                end
            else
                for n = ws.node_check
                    output{ws.idx_fnc}.histograms.(['node_',num2str(n)])=[];
                end
            end
        case 7
            %% Signal Spectra           
            output{ws.idx_fnc}.fcn_type = 'Signal Spectra';
            [ output,ws ] = search_done_fcn( 'Frame Filtering', output,ws, 0);
            if isnan(ws.sequence_num)
                [Frame, output,ws] = get_field('Frame', 'select and load .data file' ,output,ws);
            else
                [ idx ] = idx_done_fcn( 'divido nei vari nodi', output );
                % Frame Filtering
                Frame = output{idx}.Frame_node.node_1;
            end
            
            
            if ws.flag.spectra == 1
                Frame = Frame + rand(size(Frame)) - 0.5;
                ChNum = ws.ChNum;
                if (ws.flag.Tile == 1)
                    bins = linspace(0,9*max(max(Frame)),2^12); %per singolo cristallo su tile        
                else
                    bins = linspace(0,ChNum*max(max(Frame)),0.01*ChNum*2^12); %per preclinico o clinico
                end

                %plot spectra
                output{ws.idx_fnc}.spectrum = PlotSectra(Frame,ws.num_nodes,ws.Node,bins,ws.FName);
                drawnow
            else
               [node_check, output,ws] = get_field('node_check', 'select and load .data file' ,output,ws);
               for n = node_check
                   output{ws.idx_fnc}(n).spectrum = [];
               end
            end
            
        case 8
            %% Gaussian Fitting e Risoluzione Energetica
            output{ws.idx_fnc}.fcn_type = 'Gaussian Fitting e Risoluzione Energetica';
            
            [ output,ws ] = search_done_fcn( 'Signal Spectra', output,ws, 1);
            [ idx_spectra ] = idx_done_fcn( 'Signal Spectra', output );
            [ idx ] = idx_done_fcn( 'Gaussian Fitting e Risoluzione Energetica', output );

            if (ws.flag.RS == 1)
                for n = ws.node_check
                    while( menu('Choise the Spectra figure, zoom and press the "ok" button below, then select the photopeak apex', 'ok') ~= 1)
                    end
                    [xpoint, ypoint] = ginput(1);
                    fitting = monoG_v2(output{idx_spectra}.spectrum.(['node_',num2str(n)]).x, output{idx_spectra}.spectrum.(['node_',num2str(n)]).y, xpoint, ypoint);
                    output{idx}.GFitRE.(['node_',num2str(n)]).fitting = fitting;
                    output{idx}.GFitRE.(['node_',num2str(n)]).EnergyRisolution = 2.3548*fitting.c1/sqrt(2)/fitting.b1;
                    
                    title(['Energy Resolution = ', num2str(100 * round(output{idx}.GFitRE.(['node_',num2str(n)]).EnergyRisolution,4)), ' %'])
                end
            end
        case 9
            %% Images (Modified Centroid Method Reconstruction)
            output{ws.idx_fnc}.fcn_type = 'Images (Modified Centroid Method Reconstruction)';
            [Frame, output,ws] = get_field('Frame', 'select and load .data file' ,output,ws);
            [ output,ws ] = search_done_fcn( 'Signal Spectra', output,ws, 1);
            [ idx ] = idx_done_fcn( 'Images (Modified Centroid Method Reconstruction)', output );
            
            if ws.flag.images == 1
                for n=1:ws.num_nodes
                    %% seleziono la finestra energetica corrispondente al fotopicco
                    while( menu('Choise the Spectra figure, zoom and press the "ok" button below, then select the photopeak energy window', 'ok') ~= 1)
                    end
                    [q,~]=ginput(2);
                    Filt.E_min{n}=round(min(q));%channels
                    Filt.E_max{n}=round(max(q));%channels
                    if (ws.flag.multipleEW==0)
                        Filt.E_min=cell2mat(Filt.E_min);%channels
                        Filt.E_max=cell2mat(Filt.E_max);%channels
                        break
                    end
                end

                % Baseline subtraction for modified centroid method
                % Baseline value is arbitrary and has to be defined by observing the result
                % obtained in term of reconstructed image employing the Centroid Method
                for i=ws.flag.baseline
                    Filt.baseline=i;
                    output{idx}.CentroidRecostruction=PlotImages(Frame,ws.modality,Filt,ws.num_nodes,ws.Node,ws.FName);
                    drawnow
                    set(gcf, 'Name', [get(gcf, 'Name'), '___Baseline = ', num2str(i)])
                end
            %     Filt.baseline= 1.5 * mean(reshape(Frame,1,[]));%channels
            %     % plot images
            %     [output]=PlotImages(Frame,modality,Filt,num_nodes,Node,output,FName);
            %     drawnow
            end
        case 10
            %% Spettri Locali
            if (ws.flag.local_spectrum == 1)
                output{ws.idx_fnc}.fcn_type = 'Spettri Locali';
                
                [ output,ws ] = search_done_fcn( 'Images (Modified Centroid Method Reconstruction)', output,ws, 1);
                for n = ws.node_check
                    [Frame_node, output,ws] = get_field( ['Frame_node.node_',num2str(n)], 'divido nei vari nodi', output,ws );
                    binCh = linspace(0,round(max(max(Frame_node))),2^12);
                    bins = linspace(0,ws.ChNum*max(max(Frame_node)),0.01*ws.ChNum*2^12);
                    [modality, output,ws] = get_field( 'modality', 'select and load .data file', output,ws );
                    [ idx ] = idx_done_fcn( 'Images (Modified Centroid Method Reconstruction)', output );
                    
                    while( menu('Choise the Image figure, press the "ok" button below, then select with 2 input the local area', 'ok') ~= 1)
                    end
                    [p,q] = ginput(2); % prese in pixel
                    p = sort(p);
                    q = sort(q);
                    x_rec = output{idx}.CentroidRecostruction.(['node_',num2str(n)]).x_rec;
                    y_rec = output{idx}.CentroidRecostruction.(['node_',num2str(n)]).y_rec;
                    En_filter = output{idx}.CentroidRecostruction.(['node_',num2str(n)]).En_filter;
                    
                    c_bool = En_filter' & (x_rec>min(p)) & (x_rec<max(p)) & (y_rec>min(q)) & (y_rec<max(q));
                    filtri = En_filter' .* c_bool;
                    filtri(filtri==0) = [];
                    
                    subFrame = Frame_node(filtri,:);
                    [y,x] = hist( sum(subFrame,2), 75 );
                    figure, plot(x,y)                           % ssum degli eventi nell'area selezionata
                    
                    
%                     [ mCorrezione, tuttoFrame ] = correzione_spettro_locale_2( modality,output{idx}.CentroidRecostruction.(['node_',num2str(n)]),binCh,Frame_node,ws.flag.sector  );
%                     % sector = 1 -> Tile
%                     % sector = 2 -> Pixel
%                     % sector = 3 -> Complessivo
%                     % sector = 4 -> pixel Tile
%                     FrameCorretto = tuttoFrame{1,ws.flag.sector};
%                     noZero = (sum(FrameCorretto,2)==0);
%                     FrameCorretto( noZero, : ) = [];
%                     [output]=PlotSectra(FrameCorretto,1,ones(size(FrameCorretto,1),1),bins,ws.FName);

%                     % calcola RS dopo la correzione per spettri locali
%                     output.peak_fitting{1} = monoG(output.spectrum.x, output.spectrum.y, 10);
%                     output.peak_fitting{2} = 2.3548*output.peak_fitting{1}.c1/sqrt(2)/output.peak_fitting{1}.b1;
%                     title(['Energy Resolution = ', num2str(100 * output.peak_fitting{2}), ' %'])
                end   
            end
            
            
            case 19
            %% Spettri Locali 2
            if (ws.flag.local_spectrum == 1)
                output{ws.idx_fnc}.fcn_type = 'Spettri Locali 2';
                
                [ output,ws ] = search_done_fcn( 'Images (Modified Centroid Method Reconstruction)', output,ws, 1);
                [ idx ] = idx_done_fcn( 'Images (Modified Centroid Method Reconstruction)', output );
                
                for n = ws.node_check
                    [Frame_node, output,ws] = get_field( ['Frame_node.node_',num2str(n)], 'divido nei vari nodi', output,ws );
 
                    rep = 1;
                    while(rep ~= 0)
                        if (menu('Choise the Image figure, press the "ok" button below, then select with 2 input the local area', 'ok','stop') == 2)
                            rep = 0;
                        else
                            [p,q] = ginput(2); % prese in pixel
                            p = sort(p);
                            q = sort(q);
                            x_rec = output{idx}.CentroidRecostruction.(['node_',num2str(n)]).x_rec;
                            y_rec = output{idx}.CentroidRecostruction.(['node_',num2str(n)]).y_rec;
                            En_filter = output{idx}.CentroidRecostruction.(['node_',num2str(n)]).En_filter;

                            c_bool = En_filter' & (x_rec>min(p)) & (x_rec<max(p)) & (y_rec>min(q)) & (y_rec<max(q));
                            filtri = En_filter' .* c_bool;
                            filtri(filtri==0) = [];

                            subFrame = Frame_node(filtri,:);
                            binnante = @(x) (log(x) / log(10))*2^6;
                            [y,x] = hist( sum(subFrame,2), round( binnante(size(subFrame,1))) );
                            figure, plot(x,y)                           % ssum degli eventi nell'area selezionata
                            a = round([p,q]);
                            title(['from (', num2str(a(1,:)),') to (',num2str(a(2,:)),')'])
                        end
                    end
                end   
            end
        case 36
            %% Spettri Locali Spots
            output{ws.idx_fnc}.fcn_type = 'Local Spectra by Spots';
            spots = ws.spots_ROI(1:end); % ------------------------------------------------------------> needed Offset Subtraction
            [ answers ] = inputdlg({'Enter Photopeak Energy:','Enter number of bins'},'Input',1,{'122','100'});
            photopeackEnergy = str2double(answers{1});
            bins = str2double(answers{2});
            subFrame_originale = cell(numel(spots),1);
            subFrame_corretto = cell(numel(spots),1);
            for i=1:numel(spots)
%                 while( menu('Zoom and press the "ok" button below, then select the photopeak apex', 'ok') ~= 1)
%                 end
                [y,x] = hist(spots{i}.ssum, bins);
                figure('Name',['Spot number ',num2str(i),' /',num2str(numel(spots))]), plot(x,y)
                [xpoint, ypoint] = ginput(1);
                fitting = monoG_v2(x, y, xpoint, ypoint);
                spots{i}.slope = photopeackEnergy / fitting.b1;
                subFrame_originale{i,1} = spots{i}.subFrame;
                spots{i}.subFrame = spots{i}.subFrame .* spots{i}.slope;
                close gcf
                subFrame_corretto{i,1} = spots{i}.subFrame;
            end
            subFrame_originale = cell2mat(subFrame_originale);
            subFrame_corretto = cell2mat(subFrame_corretto);
            [y_originale,x_originale] = hist(sum(subFrame_originale.*spots{1}.slope,2), linspace(0,200,bins));
            [y_corretto,x_corretto] = hist(sum(subFrame_corretto,2), linspace(0,200,bins));
            h = figure();
            hold on
            plot(x_originale,y_originale,'blue')
            plot(x_corretto,y_corretto,'red')
            legend('originale', 'corretto')
            output{ws.idx_fnc}.subFrame_originale = subFrame_originale;
            output{ws.idx_fnc}.subFrame_corretto = subFrame_corretto;
            output{ws.idx_fnc}.spots = spots;
            
            % calcolo risoluzione energetica
            axesObjs = get(h, 'Children');  %axes handles
            dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
            %objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
            xdata = get(dataObjs{2}, 'XData');  %data from low-level grahics objects
            ydata = get(dataObjs{2}, 'YData');
            zdata = get(dataObjs{2}, 'ZData');

            figure, plot(xdata{1, 1}, ydata{1, 1})
            fitting = monoG(xdata{1, 1}, ydata{1, 1}, 0);
            EnergyRisolution = 2.3548*fitting.c1/sqrt(2)/fitting.b1;
            title(['risoluzione energetica = ', num2str(round(100*EnergyRisolution,2)), ' %'])
        case 11
            %% CALIBRAZIONE ENERGETICA
            if (flag.CE==1)
                FilterSpec = '*.mat';
                [filename_Ba133, pathname_Ba133] = uigetfile([DefaultPathnames.Ba133Pathname,'\',FilterSpec], 'select .mat Ba133 elaboration file', 'MultiSelect', 'on');
                DefaultPathnames.Ba133Pathname = pathname_Ba133;

                [ coba, baro ] = calibrazione_energetica( output, pathname_Ba133, filename_Ba133 );
                output = coba;
                save(['Outputs\', FName], 'output');
            end
        case 44
            %% ricavo rette di calibrazione
            output{ws.idx_fnc}.fcn_type = 'Get Calibration Lines';
            [ output,ws ] = search_done_fcn( 'select and load .data file', output,ws, 1);
            [ idx ] = idx_done_fcn( 'select and load .data file', output );
            [ idx_fnc ] = idx_done_fcn( 'Get Calibration Lines', output );

            %% set CalibrationSettings
            MinPeakDistance = 20;
            v_imp_calib = [200 300 400 500 600 700 800 900];
            initial = 100;
            th = 300;
            %
            CalibrationSettings.MinPeakDistance = MinPeakDistance;
            CalibrationSettings.v_imp_calib = v_imp_calib;
            CalibrationSettings.initial = initial;
            CalibrationSettings.th = th;
            
            %% execute the calibration
            impulsazione = 200:100:900;
            flag_exit = 0;
            while (flag_exit == 0)
                switch menu('execute:', 'execute calibration', 'add a calibration dataset', 'show Frame Istogram and Spectra', 'change values of impulsation', 'exit')
                    case 1
                        switch menu('load new dataset?', 'yes', 'no')
                            case 1
                                unisci_dati_per_calibrazione
                            case 2
                                [ idxF ] = idx_done_fcn( 'select and load .data file', output );
                                F = output{idxF}.Frame;
                        end
                        [ posizioni, coefficienti ] = new_equalizzazioneCanali( impulsazione, F );
                        output{idx_fnc}.posizioni = posizioni;
                        output{idx_fnc}.coefficienti = coefficienti;
                        %[ ~, m_corr, q_corr ] = equalizz_ch( CalibrationSettings, output{idx}.Frame, ws.ChNum, 1 );
                        uisave('coefficienti','calibrationCoefficient.mat')
                        flag_exit = 1;
                    case 2
                        [filename, pathname] = uigetfile([ws.DefaultPathnames.DataPathname,'\','*data'], 'select .data acquisition file', 'MultiSelect', 'off');
                        fid = fopen([pathname,filename],'r','b'); % per win
                        [Frame,Node,Time_stamp,modality] = openDataFile(fid);
                        output{idx}.Frame = [output{idx}.Frame; Frame];
                        ws.Node = [ws.Node; Node];
                        [ output,ws ] = do_fcn_type( 'divido nei vari nodi', output,ws );
                    case 3
%                         fid = fopen([output{idx}.pathname,output{idx}.filename],'r','b'); % per win
%                         [Frame,Node,Time_stamp,modality] = openDataFile(fid);
%                         output{idx}.Frame = Frame;
%                         ws.Node = Node;
                        [ output,ws ] = do_fcn_type( 'Signal Spectra', output,ws );
                        [ output,ws ] = do_fcn_type( 'Histograms', output,ws );
                    case 4
                        prompt = {'Enter impulsation vector:'};
                        dlg_title = 'Input';
                        num_lines = 1;
                        defaultans = {'200:100:900'};
                        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
                        impulsazione = eval(answer{1});
                    case 5
                        flag_exit = 1;
                        close all
                end
            end
            
        case 12
            %% SAVE AS .mat FILE
            save(FName, 'Frame')
        case 13
            %% Save Output
            save(['Outputs\', FName], 'output');
    end
end