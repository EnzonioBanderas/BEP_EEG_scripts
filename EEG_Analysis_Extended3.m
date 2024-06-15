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

Settings.Channels={'right-S';'right-M';'left-M';'left-S'};
Settings.State={'AA';'NREM';'REM';'QA'};
Settings.AccelerometerDimensions={'x';'y';'z'};
Settings.Bands={'delta';'theta';'alpha';'beta';'mu';'gamma'};
Settings.BandRanges=[.5,5;6,9;10,15;16,22.75;23,31.75;35,45;.5,31.75;.5,100]; %Fenzl et al. 2007, Kreuzer et al. 2015, \delta \theta \alpha \mu/\eta \beta (\gamma added)
% Settings.BandRanges=[1.5,6;6,10;10.5,15;22,30;35,45]; %Louis et al. 2004, \delta \theta \alpha \beta \gamma

% Some constants from Settings
Settings.nChannel=length(Settings.Channels);
Settings.nState=length(Settings.State);
Settings.nDim=length(Settings.AccelerometerDimensions);
Settings.nBand=size(Settings.BandRanges,1);

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
Settings.SaveFolder='C:\Users\w10-A08CFDC2D8F5\Documents\Neurologger\Analyzed\';
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
load(E_Name{E_index},'Acceleration','Electrical','header')

%% Create time vectors and general constants
Electrical.n=size(Electrical.CH1234,1);
Acceleration.n=size(Acceleration.XYZ,1);
Electrical.t=(0:1/Electrical.fs:(Electrical.n-1)/Electrical.fs)'; % s
Acceleration.t=(0:1/Acceleration.fs:(Acceleration.n-1)/Acceleration.fs)'; % s
% Some constants
Settings.nDim=size(Acceleration.XYZ,2);
Acceleration.nEpochLength=round(Settings.window_epoch*Acceleration.fs); % around Settings.window_epoch seconds
Acceleration.nEpochOverlap=round(Acceleration.nEpochLength*Settings.overlap_epoch); % around Settings.overlap_epoch*Settings.window_epoch seconds
EleAccMultiple=round(Electrical.fs/Acceleration.fs);
Electrical.nEpochLength=Acceleration.nEpochLength*EleAccMultiple;
Electrical.nEpochOverlap=Acceleration.nEpochOverlap*EleAccMultiple;

%% Calculate different parameters for each epoch
%% Compute Acceleration parameter
% Settings.nDim=1;
% Approximate dynamic acceleration by applying a high-pass filter and calculate length of vector
% approximate XYZ dynamic acceleration by applying a high-pass filter
Acceleration.dyn=[HighPassfilter(Settings.order_acc,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,1)),...
                  HighPassfilter(Settings.order_acc,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,2)),...
                  HighPassfilter(Settings.order_acc,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,3))];           
% calculate length of XYZ dynamic acceleration vector
% Acceleration.dyn=vecnorm(Acceleration.dyn,2,2);
Acceleration.dyn=sqrt(sum(Acceleration.dyn.^2,2));
% divide dynamic acceleration vector length time series into overlapping bins (1 min bins)
Acceleration.Epoch=buffer_dis(Acceleration.dyn,Acceleration.nEpochLength,Acceleration.nEpochOverlap,Acceleration.points);
Epoch.n=size(Acceleration.Epoch,2);
% calculate mean of length as parameter of bin
Epoch.DynAcc(:,1)=mean(Acceleration.Epoch,1);

Acceleration=rmfield(Acceleration,'Epoch');
%% Compute BandPower parameters (Settings.nChannel*Settings.nBand dimensions)
% Settings.nDim=Settings.nChannel*Settings.nBand

% FFT settings
FFT.nWindowLength = round(Settings.window_FFT*Electrical.fs);
FFT.nOverlap      = round(FFT.nWindowLength*Settings.overlap_FFT);
FFT.n             = FFT.nWindowLength;
FFT.BinFreq       = (0:Electrical.fs/FFT.n:Electrical.fs/2)'; % Center frequencies of DFT bins
FFT.nBinFreq      = length(FFT.BinFreq);

% preallocation
Epoch.PSD = cell(Settings.nChannel,1);
Epoch.BandPowers = cell(Settings.nChannel,1);
for i=1:Settings.nChannel % Loop through channels
Epoch.PSD{i} = zeros(Epoch.n,FFT.nBinFreq); 
Epoch.BandPowers{i} = zeros(Epoch.n,Settings.nBand);
end

% Compute band power parameter for each channel
for i=1:Settings.nChannel % loop through channels
% Filter data and divide into overlapping epochs
Electrical.Epoch = HighLowPassfilter(Settings.order_ele, [Settings.HP_ele,Settings.LP_ele], Electrical.fs, Electrical.CH1234(:,i));
Electrical.Epoch = buffer_dis(Electrical.Epoch,Electrical.nEpochLength,Electrical.nEpochOverlap,Electrical.points);
for ii=1:Settings.nBand % loop through bands
for iii=1:Epoch.n % loop through bins
    % Compute PSD and BandPowers
    [~,~,~,PSD]=spectrogram(Electrical.Epoch(:,iii),FFT.nWindowLength,FFT.nOverlap,FFT.n,Electrical.fs);
    Epoch.PSD{i}(iii,:)=mean(PSD,2);
    Epoch.BandPowers{i}(iii,ii)=bandpower2(Epoch.PSD{i}(iii,:)',FFT.BinFreq,Settings.BandRanges(ii,:));
end
end
end

%% Compute RMS, DELTA, THETA, CC and lag/delay parameters for each Epoch (and each Band and Channel)
% Settings.nDim=Settings.nChannel*Settings.nBand+2*Settings.nChannel+((Settings.nChannel^2-Settings.nChannel)/2)*Settings.nBand;
% CC settings
CC.nMaxLag       = round(Settings.maxlag_CC*Electrical.fs);
CC.nWindowLength = round(Settings.window_CC*Electrical.fs);
CC.nOverlap      = round(CC.nWindowLength*Settings.overlap_CC);
CC.nWindow       = floor((Electrical.nEpochLength-CC.nOverlap)/(CC.nWindowLength-CC.nOverlap));
CC.Lag           = (-CC.nMaxLag:CC.nMaxLag)/Electrical.fs;

% preallocation
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

for ii=1:Settings.nBand % loop through bands
    
    %%% Process and RMS
    for i=1:Settings.nChannel % Loop through channels
    Electrical.Epoch{i} = HighLowPassfilter(Settings.order_ele, Settings.BandRanges(ii,:), Electrical.fs, Electrical.CH1234(:,i)); % Filter
    Electrical.Epoch{i} = buffer_dis(Electrical.Epoch{i},Electrical.nEpochLength,Electrical.nEpochOverlap,Electrical.points); % Divide in overlapping segments
    Epoch.RMS{i}(:,ii) = rms(Electrical.Epoch{i}); % RMS
    end

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

%%% DELTA and THETA
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

%% Save Epoch as Multidimensional Dataset
Epoch_Dataset=[Epoch.BandPowers{:},... %1:32
               Epoch.RMS{:},... %33:64
               Epoch.DELTA{:},Epoch.THETA{:},Epoch.DELTA2{:},Epoch.THETA2{:},... %65:80
               Epoch.CC{triu(true(Settings.nChannel),1)},... %66:128
               Epoch.DynAcc]; %129  
%            Epoch.lag{triu(true(Settings.nChannel),1)},


%% Preallocation and settings for State Decision Model fitting
%%% Meaning of abbreviation:
% AA              = Active Awake
% (NREM||REM||QA) = Sleep
% NREM            = Non-Rapid Eye Movement
% REM             = Rapid Eye Movement
% QA              = Quiet Awake

%%% preallocation
% Probabilities and decision tree models
Epoch.p         = cell(Settings.nChannel,1);
Epoch.p_NREM_REMQA = cell(Settings.nChannel,1); % decision 2 probabilities
GMModel_NREM_REMQA=cell(Settings.nChannel,1); % decision 2 model
Epoch.ind_NREM_REMQA=cell(Settings.nChannel,1); % decision 2 index
Epoch.p_REM_QA   = cell(Settings.nChannel,1); % decision 3 probabilities
GMModel_REM_QA=cell(Settings.nChannel,1); % decision 3 model
Epoch.ind_REM_QA=cell(Settings.nChannel,1); % decision 3 index
Epoch.C         = cell(Settings.nChannel,1); % Epoch colors change with Epoch state probabilities
for i=1:Settings.nChannel
Epoch.p{i}=zeros(Epoch.n,Settings.nState);
Epoch.C{i}=zeros(Epoch.n,3);
end
% Thresholds
thresh_AA_Sleep=zeros(Settings.nChannel,1);
thresh_NREM_REMQA=zeros(Settings.nChannel,1);
thresh_REM_QA=zeros(Settings.nChannel,1);
% State logicals
Epoch.AA      = cell(Settings.nChannel,1);
Epoch.Sleep   = cell(Settings.nChannel,1);
Epoch.NREM    = cell(Settings.nChannel,1);
Epoch.REMQA   = cell(Settings.nChannel,1);
Epoch.REM     = cell(Settings.nChannel,1);
Epoch.QA      = cell(Settings.nChannel,1);
Epoch.nAA    = zeros(Settings.nChannel,1);
Epoch.nSleep   = zeros(Settings.nChannel,1);
Epoch.nNREM    = zeros(Settings.nChannel,1);
Epoch.nREMQA   = zeros(Settings.nChannel,1);
Epoch.nREM     = zeros(Settings.nChannel,1);
Epoch.nQA      = zeros(Settings.nChannel,1);
for i=1:Settings.nChannel
Epoch.AA{i}    = false(Epoch.n,1);
Epoch.Sleep{i}   = false(Epoch.n,1);
Epoch.NREM{i}    = false(Epoch.n,1);
Epoch.REMQA{i}   = false(Epoch.n,1);
Epoch.REM{i}     = false(Epoch.n,1);
Epoch.QA{i}      = false(Epoch.n,1);
end
% Plot settings
Epoch.S=Settings.S*ones(Epoch.n,1);

%% Uknown to AA/(NREM||REM||QA)
% Fit mixed Gaussian model (n=2) and calculate probabilities
GMModel_AA_Sleep = fitgmdist(Epoch.DynAcc,2,'Options',Settings.Options,'Replicates',Settings.Replicates);
[~,~,Epoch.p_AA_Sleep,~,~]=cluster(GMModel_AA_Sleep,Epoch.DynAcc);

% Order states so that AA has column 1 and Sleep has column 2
[~,ind.AA_Sleep]=sort(GMModel_AA_Sleep.mu,'descend');
for i=1:Settings.nChannel
Epoch.p{i}(:,[1,2])=Epoch.p_AA_Sleep(:,[ind.AA_Sleep(1),ind.AA_Sleep(2)]);
end

% Take the maximum probability as the state probability
[~,Epoch.ind_AA_Sleep]=max(Epoch.p{1},[],2);
for i=1:Settings.nChannel
Epoch.AA{i}=Epoch.ind_AA_Sleep==1;
Epoch.nAA(i)=sum(Epoch.AA{i});
Epoch.Sleep{i}=Epoch.ind_AA_Sleep==2;
Epoch.nSleep(i)=sum(Epoch.Sleep{i});
end

% Plot settings for each epoch
for i=1:Settings.nChannel
for iii=1:Epoch.n
    Epoch.C{i}(iii,:)=Epoch.p{i}(iii,:)*Settings.c(1:Settings.nState,:);
end
end

for i=1:Settings.nChannel

figure

subplot(1,2,1)
hold on
grid on
State={Epoch.AA;Epoch.Sleep};
nState=length(State);
Dec1=gobjects(nState,1);
for iiii=1:nState
Dec1(iiii)=scatter3(Epoch.DynAcc(State{iiii}{i}),Epoch.DELTA{i}(State{iiii}{i}),Epoch.THETA{i}(State{iiii}{i}),Epoch.S(State{iiii}{i}),Epoch.C{i}(State{iiii}{i},:),...
    'Marker',Settings.Marker{iiii},'LineWidth',.5);
end
title(Settings.Channels{i},'Interpreter','latex')
xlabel('Length Acceleration (g)','Interpreter','latex')
ylabel('NREM ratio ($\Delta=\frac{\delta*\alpha}{\beta*\gamma}$)','Interpreter','latex')
zlabel('REM ratio ($\Theta=\frac{\theta^2}{\delta*\alpha}$)','Interpreter','latex')
legend(Dec1,{'AA';'NREM'})

% plot histogram
subplot(1,2,2)
hold on
grid on
DynAcc_max=.5;
histogram(Epoch.DynAcc,'Normalization','pdf','NumBins',Settings.nBin)
% DynAcc_max=max(Epoch.DynAcc);
DynAcc=(linspace(0,DynAcc_max,Settings.N))';
plot(DynAcc,pdf(GMModel_AA_Sleep,DynAcc),'LineWidth',2)
xlim([0,DynAcc_max])
% plot threshold
gmd=GMModel_AA_Sleep;
pol(1)=-1/(2*gmd.Sigma(1))+1/(2*gmd.Sigma(2));
pol(2)=gmd.mu(1)/gmd.Sigma(1)-gmd.mu(2)/gmd.Sigma(2);
pol(3)=-gmd.mu(1)^2/(2*gmd.Sigma(1))+gmd.mu(2)^2/(2*gmd.Sigma(2))+...
log((gmd.ComponentProportion(1)*sqrt(gmd.Sigma(2)))/(gmd.ComponentProportion(2)*sqrt(gmd.Sigma(1))));
thresh_AA_Sleep(i)=max(roots(pol));
plot(ones(1,2)*thresh_AA_Sleep(i),ones(1,2).*get(gca,'YLim'),'LineWidth',2)

title('Decision 1','Interpreter','latex')
xlabel('Length Acceleration (g)','Interpreter','latex'),ylabel('Epoch count ($\%$)','Interpreter','latex')


end

%% Sleep=(NREM||REM||QA) to NREM/(REM||QA)(using NREM ratio parameter)
for i=1:Settings.nChannel
% Fit decision 2 model
DataPointsForFit=DuplicatePointsWithProbability(Epoch.DELTA{i},Epoch.p{i}(:,2),Settings.nMaxPoint);
% DataPointsForFit=Epoch.DELTA2{i};
GMModel_NREM_REMQA{i} = fitgmdist(DataPointsForFit,2,'Options',Settings.Options,'Replicates',Settings.Replicates);
[~,~,Epoch.p_NREM_REMQA{i},~,~]=cluster(GMModel_NREM_REMQA{i},Epoch.DELTA{i});

% Order
[~,ind.NREM_REMQA]=sort(GMModel_NREM_REMQA{i}.mu,'descend');
Epoch.p{i}(:,[2,3])=repmat(Epoch.p{i}(:,2),[1,2]).*Epoch.p_NREM_REMQA{i}(:,[ind.NREM_REMQA(1),ind.NREM_REMQA(2)]);
% Epoch.p{i}=Epoch.p{i}./sum(Epoch.p{i},2);

% Maximum probability
[~,Epoch.ind_NREM_REMQA{i}]=max(Epoch.p{i},[],2);
Epoch.AA{i}=Epoch.ind_NREM_REMQA{i}==1;
Epoch.nAA(i)=sum(Epoch.AA{i});
Epoch.NREM{i}=Epoch.ind_NREM_REMQA{i}==2;
Epoch.nNREM(i)=sum(Epoch.NREM{i});
Epoch.REMQA{i}=Epoch.ind_NREM_REMQA{i}==3;
Epoch.nREMQA(i)=sum(Epoch.REMQA{i});

% Plot settings for each epoch
for iii=1:Epoch.n
    Epoch.C{i}(iii,:)=Epoch.p{i}(iii,:)*Settings.c(1:Settings.nState,:);
end

figure

subplot(1,2,1)
hold on
grid on
State={Epoch.AA;Epoch.NREM;Epoch.REMQA};
nState=length(State);
Dec2=gobjects(nState,1);
for iiii=1:nState
Dec2(iiii)=scatter3(Epoch.DynAcc(State{iiii}{i}),Epoch.DELTA{i}(State{iiii}{i}),Epoch.THETA{i}(State{iiii}{i}),Epoch.S(State{iiii}{i}),Epoch.C{i}(State{iiii}{i},:),...
    'Marker',Settings.Marker{iiii},'LineWidth',1);
end
title(Settings.Channels{i},'Interpreter','latex')
xlabel('Length Acceleration (g)','Interpreter','latex')
ylabel('NREM ratio ($\Delta=\frac{\delta*\alpha}{\beta*\gamma}$)','Interpreter','latex')
zlabel('REM ratio ($\Theta=\frac{\theta^2}{\delta*\alpha}$)','Interpreter','latex')
legend(Dec2,{'AA';'NREM';'REMQA'})

% plot histogram
subplot(1,2,2)
hold on
grid on
histogram(DataPointsForFit,'Normalization','pdf','NumBins',Settings.nBin)
% histogram(Epoch.DELTA2{i}(Epoch.Sleep{i}),'Normalization','pdf')
% histogram(Epoch.DELTA2{i},'Normalization','pdf')
DELTA_max=max(Epoch.DELTA{i}(Epoch.Sleep{i}));
DELTA=(linspace(0,DELTA_max,Settings.N))';
plot(DELTA,pdf(GMModel_NREM_REMQA{i},DELTA),'LineWidth',2)
xlim([0,DELTA_max])
% plot threshold
gmd=GMModel_NREM_REMQA{i};
pol(1)=-1/(2*gmd.Sigma(1))+1/(2*gmd.Sigma(2));
pol(2)=gmd.mu(1)/gmd.Sigma(1)-gmd.mu(2)/gmd.Sigma(2);
pol(3)=-gmd.mu(1)^2/(2*gmd.Sigma(1))+gmd.mu(2)^2/(2*gmd.Sigma(2))+...
log((gmd.ComponentProportion(1)*sqrt(gmd.Sigma(2)))/(gmd.ComponentProportion(2)*sqrt(gmd.Sigma(1))));
thresh_NREM_REMQA(i)=max(roots(pol));
plot(ones(1,2)*thresh_NREM_REMQA(i),ones(1,2).*get(gca,'YLim'),'LineWidth',2)
title('Decision 2','Interpreter','latex')
xlabel('NREM ratio ($\Delta=\frac{\delta*\alpha}{\beta*\gamma}$)','Interpreter','latex'),ylabel('Epoch count (\%)','Interpreter','latex')

end

%% (REM||QA) to REM/QA (using THETA parameter)
Epoch.p{1}=Epoch.p{2}; %use probabilities of right-M for right-S
Epoch.p{4}=Epoch.p{3}; %use probabilities of left-M for left-S
% Epoch.p{1}=Epoch.p{3};
clearvars DataPointsForFit
DataPointsForFit=cell(Settings.nChannel,1);
EVA_D3=cell(Settings.nChannel,1);
for i=1:Settings.nChannel
% Fit decision 3 model
DataPointsForFit{i}=DuplicatePointsWithProbability(Epoch.THETA{i},Epoch.p{i}(:,3),Settings.nMaxPoint); %duplicate points with only REMQA
EVA_D3{i}=evalclusters(Epoch.THETA{i}(Epoch.REMQA{i}),'gmdistribution','gap','Klist',1:2);
% DataPointsForFit{i}=DuplicatePointsWithProbability(Epoch.THETA{i},sum(Epoch.p{i}(:,[1,3]),2),Settings.nMaxPoint); %duplicatepoints with awake included
% DataPointsForFit=Epoch.THETA{i};
GMModel_REM_QA{i} = fitgmdist(DataPointsForFit{i},2,'Options',Settings.Options,'Replicates',Settings.Replicates);
[~,~,Epoch.p_REM_QA{i},~,~]=cluster(GMModel_REM_QA{i},Epoch.THETA{i});

% Order
[~,ind.REM_QA]=sort(GMModel_REM_QA{i}.mu,'descend');
Epoch.p{i}(:,[3,4])=repmat(Epoch.p{i}(:,3),[1,2]).*Epoch.p_REM_QA{i}(:,[ind.REM_QA(1),ind.REM_QA(2)]);
% Epoch.p{i}=Epoch.p{i}./sum(Epoch.p{i},2);

end

% Epoch.p{2}=Epoch.p{1};
% Epoch.p{3}=Epoch.p{4};
for i=1:Settings.nChannel

% Maximum probability to create state logicals
[~,Epoch.ind_REM_QA{i}]=max(Epoch.p{i},[],2);
Epoch.AA{i}     = Epoch.ind_REM_QA{i}==1;
Epoch.nAA(i)    = sum(Epoch.AA{i});
Epoch.NREM{i}   = Epoch.ind_REM_QA{i}==2;
Epoch.nNREM(i)  = sum(Epoch.NREM{i});
Epoch.REM{i}    = Epoch.ind_REM_QA{i}==3;
Epoch.nREM(i)   = sum(Epoch.REM{i});
Epoch.QA{i}     = Epoch.ind_REM_QA{i}==4;
Epoch.nQA(i)    = sum(Epoch.QA{i});

for iii=1:Epoch.n
    Epoch.C{i}(iii,:)=Epoch.p{i}(iii,:)*Settings.c(1:Settings.nState,:);
end

%%% Plot decision 3
% Plot 1D histogram decision threshold
figure

subplot(1,2,1)
hold on
grid on
State={Epoch.AA;Epoch.NREM;Epoch.REM;Epoch.QA};
nState=length(State);
Dec3=gobjects(nState,1);
for iiii=1:nState
Dec3(iiii)=scatter3(Epoch.DynAcc(State{iiii}{i}),Epoch.DELTA{i}(State{iiii}{i}),Epoch.THETA{i}(State{iiii}{i}),Epoch.S(State{iiii}{i}),Epoch.C{i}(State{iiii}{i},:),...
    'Marker',Settings.Marker{iiii},'LineWidth',1);
end
title(Settings.Channels{i},'Interpreter','latex')
xlabel('Length Acceleration (g)','Interpreter','latex')
ylabel('NREM ratio ($\Delta=\frac{\delta*\alpha}{\beta*\gamma}$)','Interpreter','latex')
zlabel('REM ratio ($\Theta=\frac{\theta^2}{\delta*\alpha}$)','Interpreter','latex')
legend(Dec3,{'AA';'NREM';'REM';'QA'})

% plot histogram
subplot(1,2,2)
hold on
grid on
histogram(DataPointsForFit{i},'Normalization','pdf')
% histogram(Epoch.THETA2{i}(Epoch.REMQA{i}),'Normalization','pdf')
% histogram(Epoch.THETA2{i},'Normalization','pdf')
THETA_max=max(Epoch.THETA{i}(Epoch.REM{i}|Epoch.QA{i}));
THETA=(linspace(0,THETA_max,Settings.N))';
plot(THETA,pdf(GMModel_REM_QA{i},THETA),'LineWidth',2)
xlim([0,THETA_max])
% plot threshold
gmd=GMModel_REM_QA{i};
pol(1)=-1/(2*gmd.Sigma(1))+1/(2*gmd.Sigma(2));
pol(2)=gmd.mu(1)/gmd.Sigma(1)-gmd.mu(2)/gmd.Sigma(2);
pol(3)=-gmd.mu(1)^2/(2*gmd.Sigma(1))+gmd.mu(2)^2/(2*gmd.Sigma(2))+...
log((gmd.ComponentProportion(1)*sqrt(gmd.Sigma(2)))/(gmd.ComponentProportion(2)*sqrt(gmd.Sigma(1))));
thresh_REM_QA(i)=max(roots(pol));
plot(ones(1,2)*thresh_REM_QA(i),ones(1,2).*get(gca,'YLim'),'LineWidth',2)
title('Decision 3','Interpreter','latex')
xlabel('REM ratio ($\Theta=\frac{\theta^2}{\delta*\alpha}$)','Interpreter','latex'),ylabel('Epoch count ($\%$)','Interpreter','latex')

end
clearvars DataPointsForFit

%% Decision_Summary of fitted models
Decision_Summary=cell(Settings.nChannel,1);
for i=1:Settings.nChannel
    
    Decision_Summary{i}=zeros(Settings.nState,4);
    
    [~,ind.AA_Sleep]=sort(GMModel_AA_Sleep.mu,'descend');
    Decision_Summary{i}(1,1)=GMModel_AA_Sleep.mu(ind.AA_Sleep(1));
    Decision_Summary{i}(1,2)=sqrt(GMModel_AA_Sleep.Sigma(ind.AA_Sleep(1)));
    Decision_Summary{i}(1,3)=GMModel_AA_Sleep.ComponentProportion(ind.AA_Sleep(1));
    Decision_Summary{i}(1,4)=thresh_AA_Sleep(i);
    
    [~,ind.NREM_REMQA]=sort(GMModel_NREM_REMQA{i}.mu,'descend');
    Decision_Summary{i}(2,1)=GMModel_NREM_REMQA{i}.mu(ind.NREM_REMQA(1)) ;
    Decision_Summary{i}(2,2)=sqrt(GMModel_NREM_REMQA{i}.Sigma(ind.NREM_REMQA(1))) ;
    Decision_Summary{i}(2,3)=GMModel_NREM_REMQA{i}.ComponentProportion(ind.NREM_REMQA(1));
    Decision_Summary{i}(2,4)=thresh_NREM_REMQA(i);
    
    [~,ind.REM_QA]=sort(GMModel_REM_QA{i}.mu,'descend');
    Decision_Summary{i}(3,1)=GMModel_REM_QA{i}.mu(ind.REM_QA(1));
    Decision_Summary{i}(3,2)=sqrt(GMModel_REM_QA{i}.Sigma(ind.REM_QA(1)));
    Decision_Summary{i}(3,3)=GMModel_REM_QA{i}.ComponentProportion(ind.REM_QA(1));
    Decision_Summary{i}(3,4)=thresh_REM_QA(i);
    Decision_Summary{i}(4,1)=GMModel_REM_QA{i}.mu(ind.REM_QA(2));
    Decision_Summary{i}(4,2)=sqrt(GMModel_REM_QA{i}.Sigma(ind.REM_QA(2)));
    Decision_Summary{i}(4,3)=GMModel_REM_QA{i}.ComponentProportion(ind.REM_QA(2));
    Decision_Summary{i}(4,4)=thresh_REM_QA(i);
    
end

%% Combine overlapping probabilities in probability state x time matrix for left-M and right-M
P=cell(Settings.nChannel,1);
for i=1:Settings.nChannel
    P{i}=zeros(Electrical.n,Settings.nState);
    for iiii=1:Settings.nState
        P{i}(:,iiii)=buffer_inv_dis(repmat(Epoch.p{i}(:,iiii)',[Electrical.nEpochLength,1]),...
            Electrical.nEpochOverlap,'mean',Electrical.points,Electrical.n,NaN);
    end
    P{i}=P{i}./repmat(sum(P{i},2),[1,Settings.nState]);
end

%% Combine probabilities of left and right hemispheres
if EVA_D3{1}.OptimalK==2&&EVA_D3{4}.OptimalK==2
    check_D3=1;
    Electrical_P=P{1}+P{4}; % combine probabilities of hempispheres by summing right-S and left-S (M fits less well overall for decision 3)
    Electrical_P=Electrical_P./repmat(sum(Electrical_P,2),[1,Settings.nState]);
    Acceleration_P=Electrical_P(1:5:end,:);
elseif EVA_D3{1}.OptimalK==2
    check_D3=2;
    Electrical_P=P{1};
    Acceleration_P=Electrical_P(1:5:end,:);
elseif EVA_D3{4}.OptimalK==2
    check_D3=3;
    Electrical_P=P{4};
    Acceleration_P=Electrical_P(1:5:end,:);
else % revert decision 3
    check_D3=4;
    P1=P{1};
    P1(:,3)=P1(:,3)+P1(:,4);
    P1(:,4)=0;
    P4=P{4};
    P4(:,3)=P4(:,3)+P4(:,4);
    P4(:,4)=0;
    Electrical_P=P1+P4; % combine probabilities of hempispheres by summing right-S and left-S (M fits less well overall for decision 3)
    Electrical_P=Electrical_P./repmat(sum(Electrical_P,2),[1,Settings.nState]);
    Acceleration_P=Electrical_P(1:5:end,:);
    clearvars P1 P4
end

% %% Use if in decision 3 either left-S, right-S or both do not have a clear mixed gaussian distribution (n=2)
% Electrical_P(:,3)=Electrical_P(:,3)+Electrical_P(:,4);
% Electrical_P(~isnan(Electrical_P(:,1)),4)=0;
% Acceleration_P(:,3)=Acceleration_P(:,3)+Acceleration_P(:,4);
% Acceleration_P(~isnan(Acceleration_P(:,1)),4)=0;

%% Alternative: Combine Epoch probabilities of hemispheres and use combine overlap afterwards
Epoch.p_AMS_leftright=Epoch.p{1}+Epoch.p{4};
Epoch.p_AMS_leftright=Epoch.p_AMS_leftright./repmat(sum(Epoch.p_AMS_leftright,2),[1,Settings.nState]);

%% Take maximum of probabilities to create state logicals
Electrical.points_FA=logical2points(~isnan(Electrical_P(:,1)));
[~,P_FA_index]=max(Electrical_P,[],2);
P_FA_index(isnan(Electrical_P(:,1)))=NaN;
Electrical.points_AA=logical2points(P_FA_index==1);
Electrical.points_NREM=logical2points(P_FA_index==2);
Electrical.points_REM=logical2points(P_FA_index==3);
Electrical.points_QA=logical2points(P_FA_index==4);
clearvars P_FA_index

% Adjust Electrical points (probably has no effect for standard Settings)
Electrical.points_FA(:,1)     =Electrical.points_FA(:,1)-1; Electrical.points_FA=round(Electrical.points_FA/EleAccMultiple)*EleAccMultiple; Electrical.points_FA(:,1)=Electrical.points_FA(:,1)+1;
Electrical.points_AA(:,1)     =Electrical.points_AA(:,1)-1; Electrical.points_AA=round(Electrical.points_AA/EleAccMultiple)*EleAccMultiple; Electrical.points_AA(:,1)=Electrical.points_AA(:,1)+1;
Electrical.points_NREM(:,1)   =Electrical.points_NREM(:,1)-1; Electrical.points_NREM=round(Electrical.points_NREM/EleAccMultiple)*EleAccMultiple; Electrical.points_NREM(:,1)=Electrical.points_NREM(:,1)+1;
Electrical.points_REM(:,1)    =Electrical.points_REM(:,1)-1; Electrical.points_REM=round(Electrical.points_REM/EleAccMultiple)*EleAccMultiple; Electrical.points_REM(:,1)=Electrical.points_REM(:,1)+1;
Electrical.points_QA(:,1)     =Electrical.points_QA(:,1)-1; Electrical.points_QA=round(Electrical.points_QA/EleAccMultiple)*EleAccMultiple; Electrical.points_QA(:,1)=Electrical.points_QA(:,1)+1;

% Convert Electrical points to Acceleration points
Acceleration.points_FA     =[Electrical.points_FA(:,1)+EleAccMultiple-1,Electrical.points_FA(:,2)]/EleAccMultiple;
Acceleration.points_AA     =[Electrical.points_AA(:,1)+EleAccMultiple-1,Electrical.points_AA(:,2)]/EleAccMultiple;
Acceleration.points_NREM   =[Electrical.points_NREM(:,1)+EleAccMultiple-1,Electrical.points_NREM(:,2)]/EleAccMultiple;
Acceleration.points_REM    =[Electrical.points_REM(:,1)+EleAccMultiple-1,Electrical.points_REM(:,2)]/EleAccMultiple;
Acceleration.points_QA     =[Electrical.points_QA(:,1)+EleAccMultiple-1,Electrical.points_QA(:,2)]/EleAccMultiple;
    
%% Process
for i=1:Settings.nChannel
    
    Electrical.CH1234(:,i) = HighLowPassfilter(Settings.order_ele, [Settings.HP_ele,Settings.LP_ele], Electrical.fs, Electrical.CH1234(:,i));

end

%% Start of plotting
if ~exist('Plot','var')
    Plot=true;
end
if Plot

%% Plot Acceleration and Electrical data
clearvars Plot_P_data
clearvars Plot_P_truncation
Plot_points=[Electrical.points(1,1)+0*round(3600*Electrical.fs),Electrical.points(1,1)+1*round(3600*Electrical.fs)];
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
    
end

%% Compute and Plot State-Channel-Genotype-Abs/Norm PSDs
temp_points={Electrical.points_AA,Electrical.points_NREM,Electrical.points_REM,Electrical.points_QA};
for iiii=1:Settings.nState
temp_points{iiii}=temp_points{iiii}(temp_points{iiii}(:,2)-temp_points{iiii}(:,1)+1>=FFT.nWindowLength,:);
end
PSD=cell(Settings.nChannel,Settings.nState);
for i=1:Settings.nChannel
    figure
    for iiii=1:Settings.nState
        if ~isempty(temp_points{iiii})
        [PSD_matrix,t,f]=spectrogram_dis(Electrical.CH1234(:,i),temp_points{iiii},...
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
        if ~isempty(temp_points{iiii})
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
t_AA=sum(Electrical.points_AA(:,2)-Electrical.points_AA(:,1)+1)/Electrical.fs;
t_NREM=sum(Electrical.points_NREM(:,2)-Electrical.points_NREM(:,1)+1)/Electrical.fs;
t_REM=sum(Electrical.points_REM(:,2)-Electrical.points_REM(:,1)+1)/Electrical.fs;
t_QA=sum(Electrical.points_QA(:,2)-Electrical.points_QA(:,1)+1)/Electrical.fs;
    
StateSummary=[t_AA,t_NREM,t_REM,t_QA];
StateSummary=(StateSummary./repmat(sum(StateSummary,2),[1,Settings.nChannel]))*100;
PossibleFields={'p_NREM_REMQA','ind_NREM_REMQA','p_REM_QA','ind_REM_QA','C','Sleep','nAA','nSleep','nNREM','nREMQA','nREM','nQA','S','p_AA_Sleep','ind_AA_Sleep','p_3Dn4','ind_p_3Dn4','AA_3Dn4','nAA_3Dn4','NREM_3Dn4','nNREM_3Dn4','REM_3Dn4','nREM_3Dn4','QA_3Dn4','nQA_3Dn4'};
PossibleFields=PossibleFields(isfield(Epoch,PossibleFields));
Epoch=rmfield(Epoch,PossibleFields);
save([Settings.SaveFolder,Name{1},'_',Genotype{1},'_',Day{1},'_Analyzed'],'DataTable','Acceleration','Electrical','Acceleration_P','Electrical_P','Decision_Summary','Epoch','PSD','Settings','StateSummary','FFT','CC','Mouse','Epoch_Dataset','check_D3','header')
         
DataTable_cell{E_index}=DataTable;

if E_index<nE
    clearvars -except E_Name nE Settings E_index DataTable_cell FFT CC Plot Mouse check_D3 header
end

end

% save datatable and settings
DataTable=cat(1,DataTable_cell{:});
save([Settings.SaveFolder,'DataTable'],'DataTable')
save([Settings.SaveFolder,'Settings'],'Settings')