function GMModel_HD = EEG_TestFit2(E_FileName,nState,varargin)
% Plot_points_ele=round([1.08e4,1.115e4]*204.8); % Input Argument 4
% round([9e3,1e4]*204.8)
% EEG_TestFit2('10425-02_wt_D1_Analyzed.mat',4,true,round([1.08e4,1.115e4]*204.8))

%% Input
%
load(E_FileName,'Epoch','Settings','Electrical','Acceleration')
Electrical.t=(0:1/Electrical.fs:(size(Electrical.CH1234,1)-1)/Electrical.fs)';
Epoch.n=length(Epoch.DynAcc);
%
if nargin<3
    Plot=true;
else
    Plot=varargin{1};
end
if Plot
    close all
end
%
if nargin<4
    HOURS=[0,1];
    Plot_points_ele=[Electrical.points(1,1)+HOURS(1)*round(3600*Electrical.fs),Electrical.points(1,1)+HOURS(2)*round(3600*Electrical.fs)];
else
    Plot_points_ele=varargin{2};
end

%% Save Epoch as Multidimensional Dataset
Epoch_Dataset=[Epoch.BandPowers{:},... %1:32
               Epoch.RMS{:},... %33:64
               Epoch.DELTA{:},Epoch.THETA{:},Epoch.DELTA2{:},Epoch.THETA2{:},... %65:80
               Epoch.CC{triu(true(Settings.nChannel),1)},... %66:128
               Epoch.DynAcc]; %129  
%            Epoch.lag{triu(true(Settings.nChannel),1)},

%% Alternative: High dimensional fit with Epoch_Dataset
GMModel_HD = fitgmdist(Epoch_Dataset,nState,'Options',Settings.Options,'Replicates',Settings.Replicates);
[~,~,p_HD,~,~]=cluster(GMModel_HD,Epoch_Dataset);

% Take the maximum probability as the state probability
State=cell(nState,1);
legendcell=cell(nState,1);
[~,ind_HD]=max(p_HD,[],2);
for iiii=1:nState
for i=1:Settings.nChannel
State{iiii}{i}=ind_HD==iiii;
end
legendcell{iiii}=['State',num2str(iiii)];
end

% Plot settings for each epoch
C=cell(Settings.nChannel,1);
for i=1:Settings.nChannel
C{i}=zeros(Epoch.n,3);
for iii=1:Epoch.n
    C{i}(iii,:)=p_HD(iii,:)*Settings.c(1:nState,:);
end
end

for i=1:Settings.nChannel

figure
hold on
grid on
Dec1=gobjects(nState,1);
for iiii=1:nState
Dec1(iiii)=scatter3(Epoch.DynAcc(State{iiii}{i}),Epoch.DELTA{i}(State{iiii}{i}),Epoch.THETA{i}(State{iiii}{i}),5,C{i}(State{iiii}{i},:),...
    'Marker',Settings.Marker{iiii},'LineWidth',.5);
end
title(Settings.Channels{i},'Interpreter','latex')
xlabel('Length Acceleration (g)','Interpreter','latex')
ylabel('NREM ratio ($\Delta=\frac{\delta*\alpha}{\beta*\gamma}$)','Interpreter','latex')
zlabel('REM ratio ($\Theta=\frac{\theta^2}{\delta*\alpha}$)','Interpreter','latex')
legend(Dec1,legendcell)

end


if Plot
%% Combine overlapping probabilities in probability state x time matrix for left-M and right-M
Electrical_P=zeros(Electrical.n,Settings.nState);
for iiii=1:nState
    Electrical_P(:,iiii)=buffer_inv_dis(repmat(p_HD(:,iiii)',[Electrical.nEpochLength,1]),...
    Electrical.nEpochOverlap,'mean',Electrical.points,Electrical.n,NaN);
end
Electrical_P=Electrical_P./repmat(sum(Electrical_P,2),[1,Settings.nState]);
Acceleration_P=Electrical_P(1:5:end,:);

mean(Electrical_P(~isnan(Electrical_P(:,1)),:))
%% Plot Acceleration and Electrical data
Mouse=E_FileName(1:end-13);
clearvars Plot_P_data
clearvars Plot_P_truncation
nPlotLength=round(Settings.window_plot*Electrical.fs);
figure('Name',Mouse)
nSegment=size(Plot_points_ele,1);
YLIM=[-500,500];
EleAccMultiple=round(Electrical.fs/Acceleration.fs);
for j=1:nSegment
    
index=Plot_points_ele(j,1):Plot_points_ele(j,2);
[Plot_t_data,Plot_t_truncation]=buffer(Electrical.t(index),nPlotLength,1,'nodelay');
Plot_t_truncation=[Plot_t_data(end);Plot_t_truncation];
    
for i=1:Settings.nChannel
subplot(5,1,i)
hold on
grid on

    [Plot_data,Plot_truncation]=buffer(Electrical.CH1234(index,i),nPlotLength,1,'nodelay');
    Plot_truncation=[Plot_data(end);Plot_truncation];
    for iiii=1:Settings.nState
        [Plot_P_data_temp,Plot_P_truncation]=buffer(Electrical_P(index,iiii),nPlotLength,1,'nodelay');
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
xlim(Plot_points_ele/Electrical.fs)
xlabel('Time (s)','Interpreter','latex'),ylabel('EEG ($\mu V$)','Interpreter','latex')
    
end
end

clearvars Plot_P_data
clearvars Plot_P_truncation
Plot_points_acc=round([Plot_points_ele(:,1)+EleAccMultiple-1,Plot_points_ele(:,2)]/EleAccMultiple);
subplot(5,1,5)
hold on
grid on
nSegment=size(Plot_points_acc,1);
for j=1:nSegment
    
    index=Plot_points_acc(j,1):Plot_points_acc(j,2);
    
    [Plot_t_data,Plot_t_truncation]=buffer(Acceleration.t(index),nPlotLength,1,'nodelay');
    Plot_t_truncation=[Plot_t_data(end);Plot_t_truncation];
    
    [Plot_data,Plot_truncation]=buffer(Acceleration.dyn(index),nPlotLength,1,'nodelay');
    Plot_truncation=[Plot_data(end);Plot_truncation];
    for iiii=1:Settings.nState
        [Plot_P_data_temp,Plot_P_truncation]=buffer(Acceleration_P(index,iiii),nPlotLength,1,'nodelay');
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
xlim(Plot_points_ele/Electrical.fs)
    
end

subplot(5,1,1)
h=gobjects(nState,1);
for iiii=1:nState
h(iiii,1)=plot(NaN,'Color',Settings.c(iiii,:));
end
legend(h,legendcell)

%% Electrical_P replace
load(E_FileName,'Electrical_P','Acceleration_P')

%% Plot Acceleration and Electrical data
Mouse=E_FileName(1:end-13);
clearvars Plot_P_data
clearvars Plot_P_truncation

nPlotLength=round(Settings.window_plot*Electrical.fs);
figure('Name',Mouse)
nSegment=size(Plot_points_ele,1);
YLIM=[-500,500];
EleAccMultiple=round(Electrical.fs/Acceleration.fs);
for j=1:nSegment
    
index=Plot_points_ele(j,1):Plot_points_ele(j,2);
[Plot_t_data,Plot_t_truncation]=buffer(Electrical.t(index),nPlotLength,1,'nodelay');
Plot_t_truncation=[Plot_t_data(end);Plot_t_truncation];
    
for i=1:Settings.nChannel
subplot(5,1,i)
hold on
grid on

    [Plot_data,Plot_truncation]=buffer(Electrical.CH1234(index,i),nPlotLength,1,'nodelay');
    Plot_truncation=[Plot_data(end);Plot_truncation];
    for iiii=1:Settings.nState
        [Plot_P_data_temp,Plot_P_truncation]=buffer(Electrical_P(index,iiii),nPlotLength,1,'nodelay');
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
xlim(Plot_points_ele/Electrical.fs)
xlabel('Time (s)','Interpreter','latex'),ylabel('EEG ($\mu V$)','Interpreter','latex')
    
end
end

clearvars Plot_P_data
clearvars Plot_P_truncation
Plot_points_acc=round([Plot_points_ele(:,1)+EleAccMultiple-1,Plot_points_ele(:,2)]/EleAccMultiple);
subplot(5,1,5)
hold on
grid on
nSegment=size(Plot_points_acc,1);
for j=1:nSegment
    
    index=Plot_points_acc(j,1):Plot_points_acc(j,2);
    
    [Plot_t_data,Plot_t_truncation]=buffer(Acceleration.t(index),nPlotLength,1,'nodelay');
    Plot_t_truncation=[Plot_t_data(end);Plot_t_truncation];
    
    [Plot_data,Plot_truncation]=buffer(Acceleration.dyn(index),nPlotLength,1,'nodelay');
    Plot_truncation=[Plot_data(end);Plot_truncation];
    for iiii=1:Settings.nState
        [Plot_P_data_temp,Plot_P_truncation]=buffer(Acceleration_P(index,iiii),nPlotLength,1,'nodelay');
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
xlim(Plot_points_ele/Electrical.fs)
    
end

subplot(5,1,1)
h(1,1)=plot(NaN,'Color',Settings.c(1,:));
h(2,1)=plot(NaN,'Color',Settings.c(2,:));
h(3,1)=plot(NaN,'Color',Settings.c(3,:));
h(4,1)=plot(NaN,'Color',Settings.c(4,:));
legend(h,{'AA';'NREM';'REM';'QA'})


end

end
