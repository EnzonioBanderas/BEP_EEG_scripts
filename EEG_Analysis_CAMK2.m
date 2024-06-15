clearvars, close all
% Settings
Settings.Channels={'right-S';'right-M';'left-M';'left-S'};
Settings.nChannel=length(Settings.Channels);
Settings.window_FFT=2; % s
Settings.overlap_FFT=0.5; % ratio
Settings.HP_ele = 0.5; % Hz
Settings.LP_ele = 100; % Hz
Settings.order_ele = 2; % order of filter
% Settings which might be used later
Settings.window_epoch=30; % s
Settings.overlap_epoch=0.5; % ratio
Settings.window_detrend=1; % s
Settings.overlap_detrend=0.5; % ratio
Settings.window_CC=2; % s
Settings.overlap_CC=.5; % ratio
Settings.maxlag_CC=.1; % s
Settings.window_plot=1; % s
Settings.HP_acc = 1; % Hz
Settings.order_acc = 2; % order of filter
Settings.Options = statset('MaxIter',1000); % options for mixed gaussian modelling
Settings.Replicates=10;
Settings.p_value=.05; % probability
Settings.normFreq_ele=50; % Hz upper limit for normalization
Settings.nMaxPoint=100;
Settings.S=10;
Settings.State={'AA';'NREM';'REMQA'};
Settings.Bands={'delta';'theta';'alpha';'beta';'mu';'gamma';'full_Kreuzer';'full'};
Settings.BandRanges=[.5,5;6,9;10,15;16,22.75;23,31.75;35,45;.5,31.75;.5,100];
Settings.nState=length(Settings.State);
Settings.nBand=size(Settings.BandRanges,1);
Settings.nCC=(Settings.nChannel^2-Settings.nChannel)/2;

% Get folder
EF = uigetdir('','Select Experiment Folder (EF)');
cd(EF)
E=dir('*_*_*.mat');
nE=length(E);

% preallocation
Name_cell=cell(nE*Settings.nChannel,1);
Genotype_cell=cell(nE*Settings.nChannel,1);
Day_cell=cell(nE*Settings.nChannel,1);
State_cell=cell(nE*Settings.nChannel,1);
State_cell(:)={'Selected'};
Channel_cell=cell(nE*Settings.nChannel,1);
PSD_cell=cell(nE*Settings.nChannel,1);

% Start of loop through experiments
bar = waitbar(0, 'Compiling Data, please wait');
for E_index=1:nE
waitbar((E_index-1)/nE, bar, 'Compiling Data, please wait');
%% Import experiment data
% import data from .mat file
load(E(E_index).name,'Electrical','header')  

Electrical.n=size(Electrical.CH1234,1);  

% If points are not included select all
if ~isfield(Electrical,'points')
    Electrical.points=[1,Electrical.n];
end

% Temptable strings
Esplit=strsplit(E(E_index).name,'_');
E_Name=Esplit{1};
E_Genotype=Esplit{2};
E_Day=Esplit{3}(1:end-4);

% FFT settings
FFT.nWindowLength = round(Settings.window_FFT*Electrical.fs);
FFT.nOverlap      = round(FFT.nWindowLength*Settings.overlap_FFT);
FFT.n             = FFT.nWindowLength;
FFT.BinFreq       = (0:Electrical.fs/FFT.n:Electrical.fs/2)'; % Center frequencies of DFT bins
FFT.nBinFreq      = length(FFT.BinFreq);
   
%% Process
for i=1:Settings.nChannel
    Electrical.CH1234(:,i) = HighLowPassfilter(Settings.order_ele, [Settings.HP_ele,Settings.LP_ele], Electrical.fs, Electrical.CH1234(:,i));
end

%% Compute PSDs
for i=1:Settings.nChannel
    
    PSD_matrix=spectrogram_dis(Electrical.CH1234(:,i),Electrical.points,...
    FFT.nWindowLength,FFT.nOverlap,FFT.n,Electrical.fs);
    Name_cell{i+(E_index-1)*Settings.nChannel,1}=E_Name;
    Genotype_cell{i+(E_index-1)*Settings.nChannel,1}=E_Genotype;
    Day_cell{i+(E_index-1)*Settings.nChannel,1}=E_Day;
    Channel_cell{i+(E_index-1)*Settings.nChannel,1}=Settings.Channels{i};
    PSD_cell{i+(E_index-1)*Settings.nChannel,1}=mean(PSD_matrix,2);

end

end
close(bar)

DataTable=table(Name_cell,Genotype_cell,Day_cell,State_cell,Channel_cell,PSD_cell,...
                'VariableNames',{'Name','Genotype','Day','State','Channel','PSD'});
save('DataTable','DataTable')
save('Settings','Settings')