clearvars, close all

index_Kreuzer=importdata('C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Kreuzer\10425-02_wt_D1_304216-12205095_EEG1artifact_1TK_autoscored.txt');
load('C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Analyzed\G1_Automatic\10425-02_wt_D1_Analyzed.mat')

Electrical_P_cell{1}=Electrical_P;
Electrical_P_cell{1}(:,1)=Electrical_P_cell{1}(:,1)+Electrical_P_cell{1}(:,4);
Electrical_P_cell{1}(:,4)=[];
Acceleration_P_cell{1}=Acceleration_P;
Acceleration_P_cell{1}(:,1)=Acceleration_P_cell{1}(:,1)+Acceleration_P_cell{1}(:,4);
Acceleration_P_cell{1}(:,4)=[];

Plot_points=[304216,12205095];
nPoint=Plot_points(2)-Plot_points(1)+1;

nEpochLength=floor(nPoint/size(index_Kreuzer,1));
index_Kreuzer=repmat(index_Kreuzer(:,1),[1,nEpochLength]);
index_Kreuzer=index_Kreuzer(:);
index_Kreuzer_NaN=NaN(Electrical.n,1);
index_Kreuzer_NaN(304216:304216+length(index_Kreuzer)-1)=index_Kreuzer;
index_Kreuzer=index_Kreuzer_NaN;
clearvars index_Kreuzer_NaN
Electrical_P_cell{2}=NaN(Electrical.n,3);
Electrical_P_cell{2}(Plot_points(1):Plot_points(2),:)=zeros(nPoint,3);
Electrical_P_cell{2}([index_Kreuzer==1,index_Kreuzer==3,index_Kreuzer==7])=1;
Acceleration_P_cell{2}=Electrical_P_cell{2}(1:5:end,:);

%% Plot Acceleration and Electrical data

Settings.State={'Awake';'NREM';'REM'};
Settings.nState=length(Settings.State);
Settings.nChannel=4;
Settings.nState=3;
Settings.window_plot=1;
Settings.c=colormap(lines);
Settings.c([1,2],:)=Settings.c([2,1],:);
close(gcf)
Settings.Channels={'right-S';'right-M';'left-M';'left-S'};
Mouse='10425-02_wt_D1';


for TEST=1:2
    
clearvars Plot_P_data
clearvars Plot_P_truncation
HOURS=[4,5];
Plot_points=[304216+round(HOURS(1)*60*60*Electrical.fs),304216+round(HOURS(2)*60*60*Electrical.fs)];
% Plot_points=round([7000,11000]*Electrical.fs); 
% Plot_points=round([20900000,21000000]);
nPlotLength=round(Settings.window_plot*Electrical.fs);
figure('Name',Mouse)
nSegment=size(Plot_points,1);
YLIM=[-500,500];
EleAccMultiple=round(Electrical.fs/Acceleration.fs);
for j=1:nSegment
    
index=Plot_points(j,1):Plot_points(j,2);
[Plot_t_data,Plot_t_truncation]=buffer(Electrical.t(index),nPlotLength,1,'nodelay');
Plot_t_truncation=[Plot_t_data(end);Plot_t_truncation];
    
for i=1:Settings.nChannel
subplot(5,1,i)
hold on
grid on

    [Plot_data,Plot_truncation]=buffer(Electrical.CH1234(index,i),nPlotLength,1,'nodelay');
    Plot_truncation=[Plot_data(end);Plot_truncation];
    for iiii=1:Settings.nState
        [Plot_P_data_temp,Plot_P_truncation]=buffer(Electrical_P_cell{TEST}(index,iiii),nPlotLength,1,'nodelay');
        Plot_P_truncation=[Plot_P_data_temp(end);Plot_P_truncation];
        Plot_P_data(:,iiii)=[mean(Plot_P_data_temp),mean(Plot_P_truncation)];
    end
    
    nPlot=size(Plot_P_data,1);
    for jj=1:nPlot
        if jj<nPlot
            plot(Plot_t_data(:,jj),Plot_data(:,jj),'Color',Plot_P_data(jj,:)*Settings.c(1:Settings.nState,:))
        else
            plot(Plot_t_truncation,Plot_truncation,'Color',Plot_P_data(nPlot,:)*Settings.c(1:Settings.nState,:))
        end
    end
    
title(Settings.Channels{i},'Interpreter','latex')
ylim(YLIM)
xlabel('Time (s)','Interpreter','latex'),ylabel('EEG ($\mu V$)','Interpreter','latex')
    
end
end

clearvars Plot_P_data
clearvars Plot_P_truncation
Plot_points=round([Plot_points(:,1)+EleAccMultiple-1,Plot_points(:,2)]/EleAccMultiple);
subplot(5,1,5)
hold on
grid on
nSegment=size(Plot_points,1);
for j=1:nSegment
    
    index=Plot_points(j,1):Plot_points(j,2);
    
    [Plot_t_data,Plot_t_truncation]=buffer(Acceleration.t(index),nPlotLength,1,'nodelay');
    Plot_t_truncation=[Plot_t_data(end);Plot_t_truncation];
    
    [Plot_data,Plot_truncation]=buffer(Acceleration.dyn(index),nPlotLength,1,'nodelay');
    Plot_truncation=[Plot_data(end);Plot_truncation];
    for iiii=1:Settings.nState
        [Plot_P_data_temp,Plot_P_truncation]=buffer(Acceleration_P_cell{TEST}(index,iiii),nPlotLength,1,'nodelay');
        Plot_P_truncation=[Plot_P_data_temp(end);Plot_P_truncation];
        Plot_P_data(:,iiii)=[mean(Plot_P_data_temp),mean(Plot_P_truncation)];
    end
    
    nPlot=size(Plot_P_data,1);
    for jj=1:nPlot
        if jj<nPlot
            plot(Plot_t_data(:,jj),Plot_data(:,jj),'Color',Plot_P_data(jj,:)*Settings.c(1:Settings.nState,:))
        else
            plot(Plot_t_truncation,Plot_truncation,'Color',Plot_P_data(nPlot,:)*Settings.c(1:Settings.nState,:))
        end
    end
    
title('Accelerometer','Interpreter','latex')
xlabel('Time (s)','Interpreter','latex'),ylabel('Length Acceleration (g)','Interpreter','latex')
    
end




end