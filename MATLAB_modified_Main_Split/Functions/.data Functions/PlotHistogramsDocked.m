function [histograms]=PlotHistogramsDocked(modality,Frame,~,n,set_lim_dynamic,tab_group)

x_ch = 1:1:2^12; %[0:1:4095]
representation_order = reshape(reshape(1:36,6,6)',1,36);
histograms.x = zeros(size(Frame,2),numel(x_ch));
histograms.y = zeros(size(Frame,2),numel(x_ch));
if strcmp(modality,'Preclinical') 

%     figure('units','normalized','outerposition',[0 0 1 1])
%     ht = suptitle(['Node ',num2str(n),' - ASIC (1 or 2) channels']);
%     set(ht,'FontSize',18,'FontWeight','bold')
    title=strcat('Node',num2str(n));
    tab_new = uitab('Parent', tab_group, 'Title',title);
    axes('parent',tab_new)
    annotation(tab_new, 'textbox',[.0 .9 .8 .1],'String',['Node ',num2str(n),' - ASIC (1 or 2) channels'],'FitBoxToText','on')    
    for index_channel = 1:size(Frame,2)
        %subplot to represent 36 graphs for the spectra, index channel point the proper channel
        subplot(6,6,representation_order(index_channel)) 
        %eventData(:,index_channel) means: take all the rows (:) of the index_channel column in the eventData file
        %hist makes the histogram of the selected data, according to te energy axis already defined
        [y1,x1] = hist(Frame(:,index_channel),x_ch) ;
        [Hax,hLine1,hLine2]=plotyy(x1(1:set_lim_dynamic-1),y1(1:set_lim_dynamic-1),x1(set_lim_dynamic:end),y1(set_lim_dynamic:end));
        hLine1.LineWidth = 2;
        hLine1.Color = 'r';
        hLine2.Color = 'b';
        Hax(1).XAxis.Limits = [0 4250];
        Hax(2).XAxis.Limits = [0 4250];
        Hax(1).YColor= 'r';
        Hax(2).YColor= 'b';
        max_YTick = max(Hax(1).YTick);
        Hax(1).YTick=[0 max_YTick/2 max_YTick];
        Hax(1).YTickLabelRotation = 60;
%         title(sprintf('Channel %d',index_channel),'Fontsize',12)  
        text(0.6,0.9,strcat('Ch',{' '},num2str(index_channel)),'FontWeight','bold','FontSize',12,'Units','normalized')
        %save data into output
        histograms.x(index_channel,:) = x1;
        histograms.y(index_channel,:) = y1;
    end
     
elseif strcmp(modality,'Clinical') 
   
%     figure('units','normalized','outerposition',[0 0 1 1])
%     ht = suptitle(['Node ',num2str(n),' - ASIC 1 channels']);
%     set(ht,'FontSize',18,'FontWeight','bold')
    title=strcat('N',num2str(n),' -A1');
    tab_new = uitab('Parent', tab_group, 'Title',title);
    axes('parent',tab_new)
    annotation(tab_new, 'textbox',[.0 .9 .5 .1],'String',['Node ',num2str(n),' - ASIC 1 channels'],'FitBoxToText','on')
    for index_channel = 1:size(Frame,2)/2
        %subplot to represent 72 graphs for the spectra, index channel point the proper channel
        subplot(6,6,representation_order(index_channel)) 
        %eventData(:,index_channel) means: take all the rows (:) of the index_channel column in the eventData file
        %hist makes the histogram of the selected data, according to te energy axis already defined
        [y1,x1] = hist(Frame(:,index_channel),x_ch) ;
        %plot data in the center of the dynamic range
        [Hax,hLine1,hLine2]=plotyy(x1(1:set_lim_dynamic-1),y1(1:set_lim_dynamic-1),x1(set_lim_dynamic:end),y1(set_lim_dynamic:end));
        hLine1.LineWidth = 2;
        hLine1.Color = 'r';
        hLine2.Color = 'b';
        Hax(1).XAxis.Limits = [0 4250];
        Hax(2).XAxis.Limits = [0 4250];
        Hax(1).YColor= 'r';
        Hax(2).YColor= 'b';
        max_YTick = max(Hax(1).YTick);
        Hax(1).YTick=[0 max_YTick/2 max_YTick];
        Hax(1).YTickLabelRotation = 60;
%         title(sprintf('Channel %d',index_channel),'Fontsize',12)
        text(0.6,0.9,strcat('Ch',{' '},num2str(index_channel)),'FontWeight','bold','FontSize',12,'Units','normalized')
        %save data into output
        histograms.x(index_channel,:) = x1;
        histograms.y(index_channel,:) = y1;
    end
   
%      figure('units','normalized','outerposition',[0 0 1 1])
%      ht = suptitle(['Node ',num2str(n),' - ASIC 2 channels']);
%      set(ht,'FontSize',18,'FontWeight','bold')
    title=strcat('N',num2str(n),' -A2');
    tab_new = uitab('Parent', tab_group, 'Title',title);
    axes('parent',tab_new)
    annotation(tab_new, 'textbox',[.0 .9 .5 .1],'String',['Node ',num2str(n),' - ASIC 2 channels'],'FitBoxToText','on') 
    for index_channel = size(Frame,2)/2+1:size(Frame,2)
        %subplot to represent 72 graphs for the spectra, index channel point the proper channel
        subplot(6,6,representation_order(index_channel-36)) 
        %eventData(:,index_channel) means: take all the rows (:) of the index_channel column in the eventData file
        %hist makes the histogram of the selected data, according to te energy axis already defined
        [y1,x1] = hist(Frame(:,index_channel),x_ch) ;
        %plot data in the center of the dynamic range
        [Hax,hLine1,hLine2]=plotyy(x1(1:set_lim_dynamic-1),y1(1:set_lim_dynamic-1),x1(set_lim_dynamic:end),y1(set_lim_dynamic:end));
        hLine1.LineWidth = 2;
        hLine1.Color = 'r';
        hLine2.Color = 'b';
        Hax(1).XAxis.Limits = [0 4250];
        Hax(2).XAxis.Limits = [0 4250];
        Hax(1).YColor= 'r';
        Hax(2).YColor= 'b';
        max_YTick = max(Hax(1).YTick);
        Hax(1).YTick=[0 max_YTick/2 max_YTick];
        Hax(1).YTickLabelRotation = 60;
%         title(sprintf('Channel %d',index_channel),'Fontsize',12)
        text(0.6,0.9,strcat('Ch',{' '},num2str(index_channel)),'FontWeight','bold','FontSize',12,'Units','normalized')
        %save data into output
        histograms.x(index_channel,:) = x1;
        histograms.y(index_channel,:) = y1;
    end
     
end

end