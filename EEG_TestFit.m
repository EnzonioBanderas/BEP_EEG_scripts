function [GMModel,EVA] = EEG_TestFit(E_FileName,ParameterName,varargin)
% Plot_points_ele=round([1.08e4,1.115e4]*204.8); % Input Argument 4
% round([9e3,1e4]*204.8)
% EEG_TestFit('10425-02_wt_D1_Analyzed.mat','Epoch.CC{1,4}(:,2)',true,round([1.08e4,1.115e4]*204.8))
% EEG_TestFit('10425-02_wt_D1_Analyzed.mat','Epoch.BandPowers{3}(:,6)',true,round([1.08e4,1.115e4]*204.8))

if nargin<3
    Plot=true;
else
    Plot=varargin{1};
end

load(E_FileName,'Epoch','Settings','Electrical','Acceleration')

if nargin<4
    HOURS=[0,1];
    Plot_points_ele=[Electrical.points(1,1)+HOURS(1)*round(3600*Electrical.fs),Electrical.points(1,1)+HOURS(2)*round(3600*Electrical.fs)];
else
    Plot_points_ele=varargin{2};
end

% Test fit of parameter
if Plot
close all
end
Electrical.t=(0:1/Electrical.fs:(size(Electrical.CH1234,1)-1)/Electrical.fs)';
Epoch.n=length(Epoch.DynAcc);
eval(['TP=',ParameterName,';'])
% EVA=evalclusters(TP,'gmdistribution','gap','Klist',1:2);
EVA=2;
% OptimalK=EVA.OptimalK;
GMModel=fitgmdist(TP,2,'Options',Settings.Options,'Replicates',5);
[~,~,TP_pnonindexed]=cluster(GMModel,TP);
[~,INDEX]=sort(GMModel.mu,'descend');
INDEX_high=INDEX(1);
INDEX_low=INDEX(2);
TP_p=TP_pnonindexed(:,[INDEX_low,INDEX_high]);
       
COLOR  =   [0         0.4470    0.7410;   %lowcolor
            0.8500    0.3250    0.0980]; %highcolor
MARKER =   {'o';'x'};

Settings.c([1,2],:)=COLOR;

if Plot

C=zeros(Epoch.n,3);
for i=1:Epoch.n
    C(i,:)=TP_p(i,:)*COLOR;
end
for i=1:Settings.nChannel
figure
hold on
grid on
[~,index]=max(TP_p,[],2);
low=index==1;
high=index==2;
State={low,high};
nState=length(State);
Dec3=gobjects(nState,1);
for iiii=1:nState
Dec3(iiii)=scatter3(Epoch.DynAcc(State{iiii}),Epoch.DELTA{i}(State{iiii}),Epoch.THETA{i}(State{iiii}),5,C(State{iiii},:),...
'Marker',MARKER{iiii},'LineWidth',1);
end
title(Settings.Channels{i},'Interpreter','latex')
xlabel('Length Acceleration (g)','Interpreter','latex')
ylabel('NREM ratio ($\Delta=\frac{\delta*\alpha}{\beta*\gamma}$)','Interpreter','latex')
zlabel('REM ratio ($\Theta=\frac{\theta^2}{\delta*\alpha}$)','Interpreter','latex')
legend(Dec3,{['Low parameter, mean=',num2str(GMModel.mu(INDEX_low))];['High parameter, mean=',num2str(GMModel.mu(INDEX_high))]})
set(gca,'CameraPosition',[-3.0596  -66.8125    8.2449])
end

end

figure
histogram(TP,'Normalization','pdf','NumBins',.5*1e3)
hold on
grid on
XLIM=get(gca,'XLim');
TP_xpdf=linspace(XLIM(1),XLIM(2),1e5)';
TP_ypdf=pdf(GMModel,TP_xpdf);
plot(TP_xpdf,TP_ypdf,'LineWidth',2)
title([ParameterName,' test'],'Interpreter','latex')
xlabel(ParameterName,'Interpreter','latex')
ylabel('Epoch probability','Interpreter','latex')
% plot threshold
gmd=GMModel;
pol(1)=-1/(2*gmd.Sigma(1))+1/(2*gmd.Sigma(2));
pol(2)=gmd.mu(1)/gmd.Sigma(1)-gmd.mu(2)/gmd.Sigma(2);
pol(3)=-gmd.mu(1)^2/(2*gmd.Sigma(1))+gmd.mu(2)^2/(2*gmd.Sigma(2))+...
log((gmd.ComponentProportion(1)*sqrt(gmd.Sigma(2)))/(gmd.ComponentProportion(2)*sqrt(gmd.Sigma(1))));
thresh=roots(pol);
thresh=thresh(thresh>=min(gmd.mu)&thresh<=max(gmd.mu));
plot(ones(1,2)*thresh,ones(1,2).*get(gca,'YLim'),'LineWidth',2)

if Plot
%% Combine overlapping probabilities in probability state x time matrix for left-M and right-M
Electrical_P=zeros(Electrical.n,Settings.nState);
for iiii=1:2
    Electrical_P(:,iiii)=buffer_inv_dis(repmat(TP_p(:,iiii)',[Electrical.nEpochLength,1]),...
    Electrical.nEpochOverlap,'mean',Electrical.points,Electrical.n,NaN);
end
Electrical_P=Electrical_P./repmat(sum(Electrical_P,2),[1,Settings.nState]);
Acceleration_P=Electrical_P(1:5:end,:);

mean(Electrical_P(~isnan(Electrical_P(:,1)),:))

[~,ind1]=max(Electrical_P,[],2);
ind1(isnan(Electrical_P(:,1)),:)=0;
ind1=ind1==2;

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
h(1)=plot(NaN,'Color',Settings.c(1,:),'LineWidth',2);
h(2)=plot(NaN,'Color',Settings.c(2,:),'LineWidth',2);
legend(h,{['Low parameter, mean=',num2str(GMModel.mu(INDEX_low))];['High parameter, mean=',num2str(GMModel.mu(INDEX_high))]})

%% Electrical_P replace
load(E_FileName,'Electrical_P','Acceleration_P')

Settings.c([1,2],:)=Settings.c([2,1],:);

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
h=gobjects(length(Settings.State),1);
for i=1:length(Settings.State)
h(i)=plot(NaN,'Color',Settings.c(i,:),'LineWidth',2);
end
legend(h,Settings.State)
%% ind1 and ind2 agreement rate
[~,ind2]=max(Electrical_P,[],2);
ind2(isnan(Electrical_P(:,1)),:)=0;
ind2=ind2==1;

ind12=ind1(~isnan(Electrical_P(:,1)));
ind22=ind2(~isnan(Electrical_P(:,1)));

CHECK=[sum(ind12&ind22),sum(~ind12&ind22);...
       sum(ind12&~ind22),sum(~ind12&~ind22)];
CHECK=CHECK/sum(CHECK(:))
CHECK(1,1)+CHECK(2,2)

end

end

