% % Some new variables
% Settings.StateSleepString=["REM","NREM"];
% nStateSleep=length(Settings.StateSleepString);

%% OLD: 1) reshape-->buffer 2) add 2 dimensional plot delta-theta, spectrogram(time_bin)
% 4 figures per file
% Plots 1) dynamical acceleration vector length distribution with a mixed gaussian
% fit and uses relative probability thresholds to divide channel data into Awake/Sleep.
% Plots 2) dynamical acceleration vector length over time.
% Plots 3) normalized spectrograms and 
%       4) absolute PSD for Awake/Sleep for each channel
% Saves PSD data with relevant measurement information into a table to be used by Plotter.m
clearvars, close all
%% Input 
Settings.window_FFT=2; % s
Settings.overlap_FFT=0.5; % ratio
Settings.window_epoch=30; % s
Settings.overlap_epoch=0.5; % ratio
Settings.window_detrend=1; % s
Settings.overlap_detrend=0.5; % ratio
Settings.window_CC=2; % s
Settings.overlap_CC=.5; % ratio
Settings.maxlag_CC=.1; % s
Settings.window_plot=1; % s

Settings.HP_ele = 0.5; % Hz
Settings.LP_ele = 100; % Hz
Settings.order_ele = 2; % order of filter
Settings.HP_acc = 1; % Hz
Settings.order_acc = 2; % order of filter

Settings.Options = statset('MaxIter',1000); % options for mixed gaussian modelling
Settings.Replicates=10;
Settings.p_value=.05; % probability
Settings.normFreq_ele=50; % Hz upper limit for normalization
Settings.nMaxPoint=100;
Settings.S=10;

Settings.Channels={'rightS';'rightM';'leftM';'leftS'};
Settings.State={'AA';'NREM';'REMQA'};
Settings.Bands={'delta';'theta';'alpha';'beta';'mu';'gamma';'full_Kreuzer';'full'};
Settings.BandRanges=[.5,5;6,9;10,15;16,22.75;23,31.75;35,45;.5,31.75;.5,100]; %Fenzl et al. 2007, Kreuzer et al. 2015, \delta \theta \alpha \mu/\eta \beta (\gamma added)
% Settings.BandRanges=[1.5,6;6,10;10.5,15;22,30;35,45]; %Louis et al. 2004, \delta \theta \alpha \beta \gamma

% Some constants from Settings
Settings.nChannel=length(Settings.Channels);
Settings.nState=length(Settings.State);
Settings.nBand=size(Settings.BandRanges,1);
Settings.nCC=(Settings.nChannel^2-Settings.nChannel)/2;

% Select to be analyzed .mat files
[E_Name,PathName] = uigetfile('*.mat','Select the file to analyse','MultiSelect', 'on');
if ischar(E_Name) %if-statement for the case that only one file is selected
    E_Name={E_Name}; 
end
cd(PathName)
nE=length(E_Name);

if nE>1
    Plot=false;
else
    Plot=true;
end

% Folder to save analyzed data to
Settings.SaveFolder='C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Analyzed\CAMK';
% Settings.SaveFolder='C:\Users\w10-A08CFDC2D8F5\Documents\Neurologger\Analyzed\';
% Settings.SaveFolder='G:\Analyzed';

%% Plot setttings
% plot property 'Color'
Settings.c=colormap(lines); % colors to most contrasting color map
Settings.c([1,2],:)=Settings.c([2,1],:); % Switch around colors (so that AA is red)
close(gcf)
Settings.Marker={'x';'o';'p';'+';'*';'.';'s';'d';'^';'h';'v';'>';'<'};
Settings.N=1e5; % number of points for plotting pdf's
Settings.nBin=100;

% plot property 'LineStyle'
Settings.ls={'-';'--';':';'-.'};

% % All plots use LaTeX as default interpreter of text
% set(groot,'defaultTextInterpreter','latex')

% preallocate
DataTable_cell=cell(nE,1);

% Start of loop through experiments
for E_index=1:nE
Mouse=E_Name{E_index}(1:end-4);

%% Import experiment data
% import data from .mat file
load(E_Name{E_index},'Electrical','header')
   
%% Process
for i=1:Settings.nChannel
    Electrical.CH1234(:,i) = HighLowPassfilter(Settings.order_ele, [Settings.HP_ele,Settings.LP_ele], Electrical.fs, Electrical.CH1234(:,i));
end

%% Create time vectors and general constants
Electrical.n=size(Electrical.CH1234,1);
Electrical.t=(0:1/Electrical.fs:(Electrical.n-1)/Electrical.fs)'; % s
% Some constants
Electrical.nEpochLength=round(Settings.window_epoch*Electrical.fs); % around Settings.window_epoch seconds
Electrical.nEpochOverlap=round(Electrical.nEpochLength*Settings.overlap_epoch); % around Settings.overlap_epoch*Settings.window_epoch seconds
Epoch.n=sum(floor((Electrical.points(:,2)-Electrical.points(:,1)-Electrical.nEpochOverlap)/(Electrical.nEpochLength-Electrical.nEpochOverlap)));

%% Calculate different parameters for each epoch
%% Compute BandPower parameters (Settings.nChannel*Settings.nBand dimensions)
% Settings.nDim=Settings.nChannel*Settings.nBand

% FFT settings
FFT.nWindowLength = round(Settings.window_FFT*Electrical.fs);
FFT.nOverlap      = round(FFT.nWindowLength*Settings.overlap_FFT);
FFT.n             = FFT.nWindowLength;
FFT.BinFreq       = (0:Electrical.fs/FFT.n:Electrical.fs/2)'; % Center frequencies of DFT bins
FFT.nBinFreq      = length(FFT.BinFreq);

% CC settings
CC.nMaxLag       = round(Settings.maxlag_CC*Electrical.fs);
CC.nWindowLength = round(Settings.window_CC*Electrical.fs);
CC.nOverlap      = round(CC.nWindowLength*Settings.overlap_CC);
CC.nWindow       = floor((Electrical.nEpochLength-CC.nOverlap)/(CC.nWindowLength-CC.nOverlap));
CC.Lag           = (-CC.nMaxLag:CC.nMaxLag)/Electrical.fs;

% preallocation
Electrical.Epoch=cell(Settings.nChannel,1);
Epoch.PSD = cell(Settings.nChannel,1);
Epoch.BandPowers = cell(Settings.nChannel,1);
for i=1:Settings.nChannel % Loop through channels
    Epoch.PSD{i} = zeros(Epoch.n,FFT.nBinFreq); 
    Epoch.BandPowers{i} = zeros(Epoch.n,Settings.nBand);
end
Epoch.RMS=cell(Settings.nChannel,1);
Epoch.DELTA=cell(Settings.nChannel,1);
Epoch.THETA=cell(Settings.nChannel,1);
Epoch.DELTA2=cell(Settings.nChannel,1);
Epoch.THETA2=cell(Settings.nChannel,1);
Electrical.Epoch=cell(Settings.nChannel,1);
for i=1:Settings.nChannel
    Epoch.RMS{i} = zeros(Epoch.n,Settings.nBand);
end
check_matrix=triu(true(Settings.nChannel),1);
CC_dis=zeros(CC.nWindow,1);
Lag_dis=zeros(CC.nWindow,1);
Epoch.CC=cell(Settings.nChannel,Settings.nChannel);
Epoch.lag=cell(Settings.nChannel,Settings.nChannel);

% Compute band power parameter for each channel
for ii=1:Settings.nBand % loop through bands
for i=1:Settings.nChannel % loop through channels
% Filter data and divide into overlapping epochs
Electrical.Epoch{i} = HighLowPassfilter(Settings.order_ele, Settings.BandRanges(ii,:), Electrical.fs, Electrical.CH1234(:,i)); % Filter
Electrical.Epoch{i} = buffer_dis(Electrical.Epoch{i},Electrical.nEpochLength,Electrical.nEpochOverlap,Electrical.points); % Divide in overlapping segments
for iii=1:Epoch.n % loop through bins
    % Compute PSD and BandPowers
    [~,~,~,PSD]=spectrogram(Electrical.Epoch{i}(:,iii),FFT.nWindowLength,FFT.nOverlap,FFT.n,Electrical.fs);
    Epoch.PSD{i}(iii,:)=mean(PSD,2);
    Epoch.BandPowers{i}(iii,ii)=bandpower2(Epoch.PSD{i}(iii,:)',FFT.BinFreq,Settings.BandRanges(ii,:));
end
    % RMS    
    Epoch.RMS{i}(:,ii) = rms(Electrical.Epoch{i});
end

%%% Compute RMS, DELTA, THETA, CC and lag/delay parameters for each Epoch (and each Band and Channel)
% Settings.nDim=Settings.nChannel*Settings.nBand+2*Settings.nChannel+((Settings.nChannel^2-Settings.nChannel)/2)*Settings.nBand;
    %%% CC and delay
    for ROW=1:Settings.nChannel
    for COL=1:Settings.nChannel
    if check_matrix(ROW,COL)

        for iii=1:Epoch.n
            [cc]=corrgram2(Electrical.Epoch{ROW}(:,iii),Electrical.Epoch{COL}(:,iii),CC.nMaxLag,CC.nWindowLength,CC.nOverlap,Electrical.fs);
            [CC_dis,Lag_index]=max(cc);
            Lag_dis=CC.Lag(Lag_index);
            Epoch.CC{ROW,COL}(iii,ii)  = mean(CC_dis);
            Epoch.lag{ROW,COL}(iii,ii) = mean(Lag_dis);
        end

    end
    end
    end

end

%% Coherence between EEG channels
% Note that Electrical.Epoch is filtered with Settings.BandRanges(end,:)
Epoch.Coherence=cell(Settings.nChannel);
for ROW=1:Settings.nChannel
for COL=1:Settings.nChannel
if check_matrix(ROW,COL)
    Epoch.Coherence{ROW,COL} = zeros(Epoch.n,FFT.nBinFreq); 
    for iii=1:Epoch.n
        [coherence,coherenceFrequencies]=mscohere(Electrical.Epoch{ROW}(:,iii),Electrical.Epoch{COL}(:,iii),...
                                                  FFT.nWindowLength,FFT.nOverlap,FFT.n,Electrical.fs);
        Epoch.Coherence{ROW,COL}(iii,:)  = coherence;
    end
end
end
end

%% DELTA and THETA
for i=1:Settings.nChannel
    delta_alpha=Epoch.RMS{i}(:,1).*Epoch.RMS{i}(:,3);
    mu_beta=Epoch.RMS{i}(:,4).*Epoch.RMS{i}(:,5);
    theta_squared=Epoch.RMS{i}(:,2).^2;
    Epoch.DELTA{i}=delta_alpha./mu_beta;
    Epoch.THETA{i}=theta_squared./delta_alpha;
    
    delta_alpha=Epoch.BandPowers{i}(:,1).*Epoch.BandPowers{i}(:,3);
    mu_beta=Epoch.BandPowers{i}(:,4).*Epoch.BandPowers{i}(:,5);
    theta_squared=Epoch.BandPowers{i}(:,2).^2;
    Epoch.DELTA2{i}=delta_alpha./mu_beta;
    Epoch.THETA2{i}=theta_squared./delta_alpha;
end

Electrical=rmfield(Electrical,'Epoch');

%% Epochs as multidimensional dataset
% Data
Epoch_Dataset=[Epoch.BandPowers{:},... % Settings.nChannel*Settings.nBand
               Epoch.RMS{:},... % Settings.nChannel*Settings.nBand
               Epoch.DELTA{:},Epoch.THETA{:},Epoch.DELTA2{:},Epoch.THETA2{:},... % Settings.nChannel*4
               Epoch.CC{check_matrix},... % Settings.nCC*Settings.nBand
               Epoch.lag{check_matrix},... % Settings.nCC*Settings.nBand
               Epoch.PSD{:},... % Settings.nChannel*FFT.nBinFreq
               Epoch.Coherence{check_matrix}]; % Settings.nCC*FFT.nBinFreq
           
% Column names
Epoch_Dataset_BandPowers_names=cell(1,Settings.nChannel*Settings.nBand);
Epoch_Dataset_RMS_names=cell(1,Settings.nChannel*Settings.nBand);
Epoch_Dataset_Ratios_names=cell(1,Settings.nChannel*4);
Epoch_Dataset_PSD_names=cell(1,Settings.nChannel*FFT.nBinFreq);
Epoch_Dataset_CC_names=cell(1,Settings.nCC*Settings.nBand);
Epoch_Dataset_lag_names=cell(1,Settings.nCC*Settings.nBand);
Epoch_Dataset_Coherence_names=cell(1,Settings.nCC*FFT.nBinFreq);
for i=1:Settings.nChannel
    % BandPowers RMS
    for ii=1:Settings.nBand
        Epoch_Dataset_BandPowers_names{ii+(i-1)*Settings.nBand}=[Settings.Channels{i},'_',Settings.Bands{ii},'_power'];
        Epoch_Dataset_RMS_names{ii+(i-1)*Settings.nBand}=[Settings.Channels{i},'_',Settings.Bands{ii},'_RMS'];
    end
    % Ratios
    Epoch_Dataset_Ratios_names{i+0*4}=[Settings.Channels{i},'_NREMratio_RMS'];
    Epoch_Dataset_Ratios_names{i+1*4}=[Settings.Channels{i},'_REMratio_RMS'];
    Epoch_Dataset_Ratios_names{i+2*4}=[Settings.Channels{i},'_NREMratio_power'];
    Epoch_Dataset_Ratios_names{i+3*4}=[Settings.Channels{i},'_REMratio_power'];
    % PSD
    for ii=1:FFT.nBinFreq
        Epoch_Dataset_PSD_names{ii+(i-1)*FFT.nBinFreq}=[Settings.Channels{i},'_',num2str(ii),'_power'];
    end
end

count=0;
for ROW=1:Settings.nChannel
for COL=1:Settings.nChannel
if check_matrix(ROW,COL)
    count=count+1;
    % CC and lag
    for ii=1:Settings.nBand
        Epoch_Dataset_CC_names{ii+(count-1)*Settings.nBand}=[Settings.Channels{ROW},'and',Settings.Channels{COL},'_',Settings.Bands{ii},'_CC'];
        Epoch_Dataset_lag_names{ii+(count-1)*Settings.nBand}=[Settings.Channels{ROW},'and',Settings.Channels{COL},'_',Settings.Bands{ii},'_lag'];
    end
    % Coherence
    for ii=1:FFT.nBinFreq
        Epoch_Dataset_Coherence_names{ii+(count-1)*FFT.nBinFreq}=[Settings.Channels{ROW},'and',Settings.Channels{COL},'_',num2str(ii),'_coherence'];
    end
end
end
end
Epoch_Dataset_names=[Epoch_Dataset_BandPowers_names,...
                     Epoch_Dataset_RMS_names,...
                     Epoch_Dataset_Ratios_names,...
                     Epoch_Dataset_CC_names,...
                     Epoch_Dataset_lag_names,...
                     Epoch_Dataset_PSD_names,...
                     Epoch_Dataset_Coherence_names];
                 
%% Fit Gaussian mixture model
[~,Epoch_PCA,~,~,Epcoh_PCA_percentages]=pca(Epoch_Dataset);
GMModel = fitgmdist(Epoch_Dataset,Settings.nState,'Options',Settings.Options,'Replicates',Settings.Replicates,'RegularizationValue',0.1);
[~,~,Epoch.p,~,~]=cluster(GMModel,Epoch_Dataset);
statePermutation=zeros(Settings.nState,1);

iRatio=2*Settings.nChannel*Settings.nBand+1:2*Settings.nChannel*Settings.nBand+1+4*4;
[~,stateMax]=max(GMModel.mu(:,iRatio));

iNREM=[1:Settings.nChannel,2*Settings.nChannel+1:3*Settings.nChannel];
statePermutation(2)=mode(stateMax(iNREM)); % NREM

iREM=[Settings.nChannel+1:2*Settings.nChannel,3*Settings.nChannel+1:4*Settings.nChannel];
temp=stateMax(iREM);
temp=temp(temp~=statePermutation(2));
if isempty(temp)
    statePermutation(3)=mode(stateMax(stateMax~=statePermutation(2)));
else
    statePermutation(3)=mode(temp); % REM
end

iPermutation=1:Settings.nState;
statePermutation(1)=iPermutation(iPermutation~=statePermutation(2)&iPermutation~=statePermutation(3));

Epoch.p=Epoch.p(:,statePermutation);

%% Plot cluster
[~,Epoch.p_index]=max(Epoch.p,[],2);
Epoch.p_states=true(Epoch.n,Settings.nState);
for i=1:Settings.nState
    Epoch.States(:,i)=Epoch.p_index==i;
end

% Color for each epoch depends on probabilities
for iii=1:Epoch.n
    Epoch.C(iii,:)=Epoch.p(iii,:)*Settings.c(1:Settings.nState,:);
end

for i=1:Settings.nChannel
    figure
    hold on
    grid on
    plotObjects=gobjects(Settings.nState,1);
    for iiii=1:Settings.nState
        plotObjects(iiii)=scatter3(Epoch.BandPowers{i}(Epoch.States(:,iiii),6),Epoch.DELTA{i}(Epoch.States(:,iiii)),Epoch.THETA{i}(Epoch.States(:,iiii)),Settings.S,Epoch.C(Epoch.States(:,iiii),:),...
        'Marker',Settings.Marker{iiii},'LineWidth',.5);
    end
    title(Settings.Channels{i},'Interpreter','latex')
    xlabel('$\gamma$ power ($\mu V^2$)','Interpreter','latex')
    ylabel('NREM ratio ($\Delta=\frac{\delta*\alpha}{\beta*\gamma}$)','Interpreter','latex')
    zlabel('REM ratio ($\Theta=\frac{\theta^2}{\delta*\alpha}$)','Interpreter','latex')
    legend(plotObjects,Settings.State)
end

%% Save Epoch_Dataset in MATLAB table form and .xlsx
Epoch_Dataset=[Epoch_Dataset,Epoch.p,[FFT.BinFreq;zeros(Epoch.n-FFT.nBinFreq,1)]];
Epoch_Dataset_names=[Epoch_Dataset_names,Settings.State',{'CenterFrequenciesInHz'}];
Epoch_Dataset_Table=array2table(Epoch_Dataset,'VariableNames',Epoch_Dataset_names);
writetable(Epoch_Dataset_Table,[Settings.SaveFolder,filesep,Mouse,'_EpochDataset.xlsx'])

%% Combine overlapping epoch probabilities into probability vector for each point
P=zeros(Electrical.n,Settings.nState);
for iiii=1:Settings.nState
    P(:,iiii)=buffer_inv_dis(repmat(Epoch.p(:,iiii)',[Electrical.nEpochLength,1]),...
        Electrical.nEpochOverlap,'mean',Electrical.points,Electrical.n,NaN);
end
P=P./repmat(sum(P,2),[1,Settings.nState]);

%% Take maximum of probabilities to create state logicals
Electrical.points=logical2points(~isnan(P(:,1)));
[~,P_index]=max(P,[],2);
P_index(isnan(P(:,1)))=NaN;
StatePoints=cell(Settings.nState,1);
for i=1:Settings.nState
    StatePoints{i}=logical2points(P_index==i);
end
clearvars P_index

%% Start of plotting
if ~exist('Plot','var')
    Plot=true;
end
if Plot

%% Plot Electrical data
% there has to be at least one continuous hour of data
clearvars Plot_P_data
clearvars Plot_P_truncation
Plot_points=[Electrical.points(1,1)+0*round(3600*Electrical.fs),Electrical.points(1,1)+1*round(3600*Electrical.fs)];
% Plot_points=round([7000,11000]*Electrical.fs); 
% Plot_points=round([20900000,21000000]);
nPlotLength=round(Settings.window_plot*Electrical.fs);
figure('Name',Mouse)
nSegment=size(Plot_points,1);
YLIM=[-500,500];
for j=1:nSegment
    
index=Plot_points(j,1):Plot_points(j,2);
[Plot_t_data,Plot_t_truncation]=buffer(Electrical.t(index),nPlotLength,1,'nodelay');
Plot_t_truncation=[Plot_t_data(end);Plot_t_truncation];
    
for i=1:Settings.nChannel
subplot(Settings.nChannel,1,i)
hold on
grid on

    [Plot_data,Plot_truncation]=buffer(Electrical.CH1234(index,i),nPlotLength,1,'nodelay');
    Plot_truncation=[Plot_data(end);Plot_truncation];
    for iiii=1:Settings.nState
        [Plot_P_data_temp,Plot_P_truncation]=buffer(P(index,iiii),nPlotLength,1,'nodelay');
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

%% Compute and Plot State-Channel-Genotype-Abs/Norm PSDs
for iiii=1:Settings.nState
StatePoints{iiii}=StatePoints{iiii}(StatePoints{iiii}(:,2)-StatePoints{iiii}(:,1)+1>=FFT.nWindowLength,:);
end
PSD=cell(Settings.nChannel,Settings.nState);
for i=1:Settings.nChannel
    figure
    for iiii=1:Settings.nState
        if ~isempty(StatePoints{iiii})
        [PSD_matrix,t,f]=spectrogram_dis(Electrical.CH1234(:,i),StatePoints{iiii},...
         FFT.nWindowLength,FFT.nOverlap,FFT.n,Electrical.fs);
     subplot(4,1,iiii)
     imagesc(t,f,PSD_matrix)
     PSD{i,iiii}=mean(PSD_matrix,2);
     title([Mouse,' ',Settings.Channels{i},' ',Settings.State{iiii}],'Interpreter','none')
     xlabel('Time (s)','Interpreter','latex'),ylabel('Frequency (Hz)','Interpreter','latex')
     set(gca,'YDir','normal')
     set(gca,'CLim',[0,mean(max(PSD_matrix))])
     ylim([0,20])
        end
    end
end

for i=1:Settings.nChannel
count=0;
    figure
    hold on
    grid on
%     plot_h=gobjects(Settings.nChannel,1);
    for iiii=1:Settings.nState
        if ~isempty(StatePoints{iiii})
            count=count+1;
    plot_h(count)=plot(f,PSD{i,iiii},'Color',Settings.c(iiii,:),'LineWidth',2);
    legend_cell{count}=Settings.State{iiii};
        end
    end
    title([Mouse,' ',Settings.Channels{i}],'Interpreter','none')
    xlabel('Frequency (Hz)','Interpreter','latex'),ylabel('Power ($\mu V^2$)','Interpreter','latex')
    xlim([0,50])
    legend(plot_h,legend_cell)
end


end % end of plotting

%% Convert variables to DataTable 
MouseSplit=strsplit(Mouse,'_');
Name=MouseSplit(1);
Genotype=MouseSplit(2); 
Day=MouseSplit(3);

nPSD_cell=Settings.nChannel*Settings.nState;

Name_cell=cell(nPSD_cell,1);
Genotype_cell=Name_cell;
Day_cell=Name_cell;
Channel_cell=Name_cell;
State_cell=Name_cell;
PSD_table_cell=cell(nPSD_cell,1);
Electrical_fs=zeros(nPSD_cell,1);

Name_cell(:)=Name;
Genotype_cell(:)=Genotype;
Day_cell(:)=Day;

count=0;
for d=1:Settings.nChannel % Channels
for dd=1:Settings.nState % States
    
    count=count+1;    
        
    Channel_cell(count)=Settings.Channels(d);
    State_cell(count)=Settings.State(dd);
    PSD_table_cell(count)=PSD(d,dd);
    Electrical_fs(count)=Electrical.fs;

end
end

TableNames={'Name','Genotype','Day','Channel','State','Power','fs'};
DataTable=table(Name_cell,Genotype_cell,Day_cell,Channel_cell,State_cell,...
                PSD_table_cell,Electrical_fs,...
                'VariableNames',TableNames);

%% Clear and save
StateSummary=zeros(Settings.nState,2);
for i=1:Settings.nState
    StateSummary(i,1)=sum((StatePoints{i}(:,2)-StatePoints{i}(:,1)))/Electrical.fs;
end
StateSummary(:,2)=StateSummary(:,1)/sum(StateSummary(:,1));

PossibleFields={'p_NREM_REMQA','ind_NREM_REMQA','p_REM_QA','ind_REM_QA','C','Sleep','nAA','nSleep','nNREM','nREMQA','nREM','nQA','S','p_AA_Sleep','ind_AA_Sleep','p_3Dn4','ind_p_3Dn4','AA_3Dn4','nAA_3Dn4','NREM_3Dn4','nNREM_3Dn4','REM_3Dn4','nREM_3Dn4','QA_3Dn4','nQA_3Dn4'};
PossibleFields=PossibleFields(isfield(Epoch,PossibleFields));
Epoch=rmfield(Epoch,PossibleFields);
save([Settings.SaveFolder,Name{1},'_',Genotype{1},'_',Day{1},'_Analyzed'],'DataTable','Electrical','P','PSD','Settings','FFT','CC','Mouse','Epoch_Dataset_Table','header','StatePoints','StateSummary')
         
DataTable_cell{E_index}=DataTable;

if E_index<nE
    clearvars -except E_Name nE Settings E_index DataTable_cell FFT CC Plot Mouse check_D3 header
end

end

% save datatable and settings
DataTable=cat(1,DataTable_cell{:});
save([Settings.SaveFolder,'DataTable'],'DataTable')
save([Settings.SaveFolder,'Settings'],'Settings')