function [ ADC_ch_corr, m_corr, q_corr ] = equalizz_ch( CalibrationSettings, eventDatar, Phys_order, flag_calibfigure )

%% IMPORTANTE DA SETTARE
MinPeakDistance = CalibrationSettings.MinPeakDistance;
v_imp_calib = CalibrationSettings.v_imp_calib;
initial = CalibrationSettings.initial;
th = CalibrationSettings.th;

%% per distinguere Preclinico da Clinico
ChNum = numel(Phys_order);
num_imp=numel(v_imp_calib);
dinamica_ADC=2^12;
x_locs=zeros(ChNum,num_imp);
x_locs_val_medio=zeros(1,num_imp);
span=5; %span per smooth


%regression coefficients
N=ChNum+1;
r=zeros(N,1);
m=r; b=r;
m_corr=r; q_corr=r;


    x_ch = [1:1:dinamica_ADC];
    
    figure('name','Initial Hist 1-36','units','normalized','outerposition',[0 0 1 1],'numbertitle','off')
    for index_channel = 1:36
        
        %Subplot to represent 72 graphs for the spectra. Index channel point
        %the proper channel
        subplot(6,6,index_channel) 
        %eventData(:,index_channel) means: take all the rows (:) of the
        %index_channel column in the eventData file
        %Hist makes the histogram of the selected data, according to te energy
        %axis already defined
        [y1,x1] = hist(eventDatar(:,index_channel),x_ch) ;
        y1=smooth(y1,span);
        plot(x1(10:end),y1(10:end))
%        title(sprintf('Channel %d',Phys_order(index_channel)),'Fontsize',12);
         
        [pks1,locs1]=findpeaks(y1,'MinPeakDistance',MinPeakDistance);
        hold on
        %delete data under threshold
        locs1(pks1<th)=[];
        %locs(1)=[];
        pks1(pks1<th)=[];
        pks1(locs1<initial)=[]; %dato che x1 parte da 10 in su fino a 4096
        locs1(locs1<initial)=[];
        plot(x1(locs1),pks1,'rx')
        xlim([0 dinamica_ADC]);
        x_locs(index_channel,:)=x1(locs1);  
        
    end
    
    if (flag_calibfigure==0)
        close 'Initial Hist 1-36'
    end
    
    if (ChNum == 72)
        figure('name','Initial Hist 37-72','units','normalized','outerposition',[0 0 1 1],'numbertitle','off')
        for index_channel = ((N-1)/2)+1:(N-1)
            %Subplot to represent 72 graphs for the spectra. Index channel point
            %the proper channel
            subplot(6,6,index_channel-(N-1)/2)
            %eventData(:,index_channel) means: take all the rows (:) of the
            %index_channel column in the eventData file
            %Hist makes the histogram of the selected data, according to te energy
            %axis already defined
            [y2,x2] = hist(eventDatar(:,index_channel),x_ch);
            y2=smooth(y2,span);
            plot(x2(10:end),y2(10:end))
            title(sprintf('Channel %d',Phys_order(index_channel)-36),'Fontsize',12);

            [pks2,locs2]=findpeaks(y2,'MinPeakDistance',MinPeakDistance);
            %delete data under threshold
            locs2(pks2<th)=[];
            pks2(pks2<th)=[];
            pks2(locs2<initial)=[]; %dato che x2 parte da 10 in su fino a 4096
            locs2(locs2<initial)=[];
            hold on
            plot(x2(locs2),pks2,'rx')
            xlim([0 dinamica_ADC]);
            x_locs(index_channel,:)=x2(locs2);
        end

        if (flag_calibfigure==0)
            close 'Initial Hist 37-72'
        end
    end

%% plot channel curves
 figure('name','Curves','numbertitle','off')
for index_channel = 1:(N-1)
   
    hold on
     plot(v_imp_calib, x_locs(index_channel,:))
     plot(v_imp_calib, x_locs(index_channel,:), 'o')
   
end

if (flag_calibfigure==0)
    close Curves
end

%% plot channel regression curves compared with real calibration values

y=zeros(N,num_imp);
y_new=y;

figure('name','Regression Curves','numbertitle','off')
for index_channel =1:(N-1)
    hold on
%     fitting = fit(v_imp_calib, x_locs(index_channel,:),'poly1');
%     plot([0, v_imp_calib, 1000], fitting([0, v_imp_calib, 1000]));
    
    [r(index_channel,1), m(index_channel,1), b(index_channel,1)]=regression(v_imp_calib, x_locs(index_channel,:));
     y(index_channel,:)=m(index_channel,1)*v_imp_calib + b(index_channel,1)*ones(1,8);
     plot(v_imp_calib, y(index_channel,:));
     scatter(v_imp_calib, x_locs(index_channel,:))
end
if (flag_calibfigure==0)
    close 'Regression Curves'
end

%% mean values calculation and mean regression curve

for index_column=1:num_imp
    x_locs_val_medio(index_column)=sum(x_locs(:,index_column))/(N-1);
end

[r(N,1), m(N,1), b(N,1)]=regression(v_imp_calib, x_locs_val_medio);
y(N,:)=m(N,1)*v_imp_calib + b(N,1)*ones(1,8);
hold on
% plot([0 v_imp_calib 1000], y(N,:), 'black', 'LineWidth',3);
 figure('name','Retta Media','numbertitle','off')
 hold on
 plot(v_imp_calib,y(N,:))
 scatter(v_imp_calib,x_locs_val_medio)
 hold off
 if (flag_calibfigure==0)
    close 'Retta Media'
 end
 
 %% per riportare le curve regressive dei singoli canali a quella della media (equalizzazione):

figure('name','Corrected Curves','numbertitle','off')
for index_channel=1:(N-1)
    
    m_corr(index_channel,1)=m(N,1)/m(index_channel,1);
    q_corr(index_channel,1)=b(N,1)-m(N,1)/m(index_channel,1)*b(index_channel,1);
    Q=ones(1,num_imp)*q_corr(index_channel,1);
    y_new(index_channel,:)=m_corr(index_channel,1)*y(index_channel,:) + Q;
    
    hold on   
    plot(v_imp_calib, y_new(index_channel,:))
end

if (flag_calibfigure==0)
    close 'Corrected Curves'
end

%% ADC_ch_corr vs ADC_ch_orig (= x_ch):

ADC_ch_orig=(x_ch)';
ADC_ch_corr=zeros(length(ADC_ch_orig),N-1);
figure('name','original ADC chs vs correct ADC chs','numbertitle','off')
for index_channel=1:(N-1)
    ADC_ch_corr(:,index_channel)=m_corr(index_channel,1)*ADC_ch_orig + q_corr(index_channel,1)*ones(length(ADC_ch_orig),1);
    hold on
    plot(ADC_ch_orig, ADC_ch_corr(:,index_channel))
end
if (flag_calibfigure==0)
    close 'original ADC chs vs correct ADC chs'
end

%% istogrammi equalizzati:

disp('Istogrammi Equalizzati:');

% riordinando le x, quindi i valori dei canali ADC:

 figure('name','Equalized Histograms 1-36','numbertitle','off','units','normalized','outerposition',[0 0 1 1])
for index_channel = 1:36
    
    %Subplot to represent 72 graphs for the spectra. Index channel point
    %the proper channel
    subplot(6,6,index_channel)
    %eventData(:,index_channel) means: take all the rows (:) of the
    %index_channel column in the eventData file
    %Hist makes the histogram of the selected data, according to te energy
    %axis already defined
    [y_eq_1,x_ch] = hist(eventDatar(:,index_channel),x_ch) ;
    y_eq_1=smooth(y_eq_1,span);
    x_eq_1=m_corr(index_channel,1)*x_ch + q_corr(index_channel,1)*ones(1,dinamica_ADC);
    plot(x_eq_1(10:end),y_eq_1(10:end))
    title(sprintf('Channel %d',Phys_order(index_channel)),'Fontsize',12);
    
    [pks_eq_1,locs_eq_1]=findpeaks(y_eq_1,'MinPeakDistance',MinPeakDistance);
    hold on
    %delete data under threshold
    locs_eq_1(pks_eq_1<th)=[];
    pks_eq_1(pks_eq_1<th)=[];
    pks_eq_1(locs_eq_1<10)=[]; %dato che x1 parte da 10 in su fino a 4096
    locs_eq_1(locs_eq_1<10)=[];
    plot(x_eq_1(locs_eq_1),pks_eq_1,'rx')
    xlim([0 dinamica_ADC]);
    x_locs_eq(index_channel,:)=x_eq_1(locs_eq_1);
    
end

if (flag_calibfigure==0)
    close 'Equalized Histograms 1-36'
end
 
if (ChNum == 72)
     figure('name','Equalized Histograms 37-72','numbertitle','off','units','normalized','outerposition',[0 0 1 1])
     for index_channel = ((N-1)/2)+1:(N-1)

         %Subplot to represent 72 graphs for the spectra. Index channel point
         %the proper channel
         subplot(6,6,index_channel-36)
         %eventData(:,index_channel) means: take all the rows (:) of the
         %index_channel column in the eventData file
         %Hist makes the histogram of the selected data, according to te energy
         %axis already defined
         [y_eq_2,x_ch] = hist(eventDatar(:,index_channel),x_ch);
         x_eq_2=m_corr(index_channel,1)*x_ch + q_corr(index_channel,1)*ones(1,dinamica_ADC);
         y_eq_2=smooth(y_eq_2,span);
         plot(x_eq_2(10:end),y_eq_2(10:end))
         title(sprintf('Channel %d',Phys_order(index_channel)-36),'Fontsize',12);

         [pks_eq_2,locs_eq_2]=findpeaks(y_eq_2,'MinPeakDistance',MinPeakDistance);
         %delete data under threshold
         locs_eq_2(pks_eq_2<th)=[];
         pks_eq_2(pks_eq_2<th)=[];
         pks_eq_2(locs_eq_2<10)=[]; %dato che x2 parte da 10 in su fino a 4096
         locs_eq_2(locs_eq_2<10)=[];
         hold on
         plot(x_eq_2(locs_eq_2),pks_eq_2,'rx')
         xlim([0 dinamica_ADC]);
         x_locs_eq(index_channel,:)=x_eq_2(locs_eq_2);
     end
     if (flag_calibfigure==0)
        close 'Equalized Histograms 37-72'
     end
end

%% plot channel curves equalizzate
figure('name','Equalized Curves','numbertitle','off')
for index_channel = 1:(N-1)
   
    hold on
    plot(v_imp_calib, x_locs_eq(index_channel,:))
    plot(v_imp_calib, x_locs_eq(index_channel,:), 'o')
   
end

if (flag_calibfigure==0)
    close 'Equalized Curves'
end

%% calcolo dell'errore

CoefficientVariation = 100 * sqrt(var(y(1:(end-1),:))) ./ y(end,:);
figure('name','Coefficient of Variation','numbertitle','off'), hold on
plot(v_imp_calib, CoefficientVariation)
scatter(v_imp_calib, CoefficientVariation, 'filled')

legend('Coefficient of Variation [%] = 100*sd/mean')
xlabel('Impulse Values')
ylabel('Coefficient of Variation')

if (flag_calibfigure==0)
    close 'Coefficient of Variation'
end



end


