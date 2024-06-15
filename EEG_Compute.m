% Compute coherence, correlation, fft etc.
% Define settings for fitting here
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

% Check if Analyzed file was selected, if so give error
for i=1:nE
    if contains(E_Name{i},'Analyzed')
        error('*_Analyzed.mat was selected, select file without Analyzed instead.')
    end
end

% Folder to save analyzed data to
Settings.SaveFolder='C:\Users\w10-A08CFDC2D8F5\Documents\Neurologger\Analyzed\G6\';
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

%% Clear and save
PossibleFields={'p_NREM_REMQA','ind_NREM_REMQA','p_REM_QA','ind_REM_QA','C','Sleep','nAA','nSleep','nNREM','nREMQA','nREM','nQA','S','p_AA_Sleep','ind_AA_Sleep','p_3Dn4','ind_p_3Dn4','AA_3Dn4','nAA_3Dn4','NREM_3Dn4','nNREM_3Dn4','REM_3Dn4','nREM_3Dn4','QA_3Dn4','nQA_3Dn4'};
PossibleFields=PossibleFields(isfield(Epoch,PossibleFields));
Epoch=rmfield(Epoch,PossibleFields);
save([Settings.SaveFolder,Mouse,'_Analyzed'],'Acceleration','Electrical','Epoch','Settings','FFT','CC','Mouse','Epoch_Dataset','header','-v7.3')
if E_index<nE
    clearvars -except E_Name nE Settings E_index FFT CC Mouse header
end

end
save([Settings.SaveFolder,'Settings'],'Settings')