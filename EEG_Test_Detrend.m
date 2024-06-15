clearvars, close all
%% Input 
Settings.window=2; % s
Settings.overlap_ratio=0.5; % ratio
Settings.HP_ele = 1; % Hz
Settings.LP_ele = 100; % Hz
Settings.HP_acc = 1; % Hz
Settings.p_value=.03; % probability
Settings.normFreq_ele=50; % Hz upper limit for normalized spectrogram
Settings.timeBin=60; % s
Settings.timeSelect=4*60*60; % s
Settings.StateString=["Awake","Sleep"]; % states

nState=length(Settings.StateString); 

% Select to be analyzed .mat files
[FileName,PathName] = uigetfile('*.mat','Select the file to analyse','MultiSelect', 'on');
if ischar(FileName) %if-statement for the case that only one file is selected
    FileName={FileName}; 
end
cd(PathName)

%% Plot setttings
% plot property 'Color'
Settings.c=colormap(lines); % colors to most contrasting color map
close all

% Switch colors
Settings.c([1,2],:)=Settings.c([2,1],:);

% plot property 'LineStyle'
Settings.ls=["-";"--";":";"-."];

%% Import select logicals
load('DataTable_Logical') % File created by Select_Truncation.m

nExperiments=length(FileName);
DataTable_cell=cell(nExperiments,1);
% Start of loop through experiments
for experiment_index=1:nExperiments
Mouse=FileName{experiment_index}(1:end-4);

% Index for DataTable_logical which is used for the truncation
temp.select_index=find(strcmp(Mouse,DataTable_Logical.Experiment_Name));

%% Import experiment data
% import data from .mat file
load(FileName{experiment_index},'Acceleration','Electrical')

%% Constants and Truncation of data depending on EEG_Select.m output and length of bins
% create time vectors (to also truncate)
Electrical.t=0:1/Electrical.fs:(size(Electrical.CH1234,1)-1)/Electrical.fs; % s
Acceleration.t=0:1/Acceleration.fs:(size(Acceleration.XYZ,1)-1)/Acceleration.fs; % s

% truncate acceleration data with logical made by Select.m
Acceleration.XYZ=Acceleration.XYZ(DataTable_Logical.Select_Acceleration{temp.select_index},:);
Acceleration.t=Acceleration.t(DataTable_Logical.Select_Acceleration{temp.select_index});
% truncate electrical data
Electrical.CH1234=Electrical.CH1234(DataTable_Logical.Select_Electrical{temp.select_index},:);
Electrical.t=Electrical.t(DataTable_Logical.Select_Electrical{temp.select_index});

% additional small truncation with modulus so that we have whole time bins
% (e.g. 20-60 s)
% truncate acceleration data
Acceleration.nBin=round(Settings.timeBin*Acceleration.fs);
temp.mod=mod(length(Acceleration.t),Acceleration.nBin);
Acceleration.XYZ=Acceleration.XYZ(1:end-temp.mod,:);
Acceleration.t=Acceleration.t(1:end-temp.mod);
% truncate electrical data
Electrical.nBin=round(Settings.timeBin*Electrical.fs);
temp.mod=mod(length(Electrical.t),Electrical.nBin);
Electrical.CH1234=Electrical.CH1234(1:end-temp.mod,:);
Electrical.t=Electrical.t(1:end-temp.mod);

% some constants for acceleration data
Acceleration.n=length(Acceleration.t);
Acceleration.nCol=Acceleration.n/Acceleration.nBin;
Acceleration.nDim=size(Acceleration.XYZ,2);
% some constants for electrical data
Electrical.n=length(Electrical.t);
Electrical.nCol=Electrical.n/Electrical.nBin;
Electrical.nChannel=size(Electrical.CH1234,2);

end
    