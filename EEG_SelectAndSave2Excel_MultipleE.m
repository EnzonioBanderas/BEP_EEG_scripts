% Enzo Nio
% enzonio@hotmail.com
% 08-03-2019
% This script lets you select a part of the neurologger recording. The
% recording is shown as EEG trace of 4 channels and 1 plot of the dynamic
% acceleration measured by the accelerometer. After selection spectra are
% computed and written to an excel file which has the exact date and time
% in its file name. The Settings can be changed in the "Define Settings and
% Channels" section. Particularly "Settings.Plot" should be set to "true" or
% "false" depending on whether you want plots to be generated for each file
% selected.
%
% Steps to use this script:
% 1. Select .mat file converted from .edf (probably best to name it as "info1_info2_info3.mat")
% 2. In the plot of EEG and accelerometer data, click 2 times to select the time in between for analysis.
% 3. Spectra are saved in "SpectraSelectionAnalysis day-month-year hour-minute-second.xlsx".
%
% Note: It might take very long
clearvars, close all
%% Define Settings and Channels 
Settings.window=2; % s
Settings.overlap=0.5; % ratio
Settings.HP_ele = 0.1; % Hz
Settings.LP_ele = 100; % Hz
Settings.order_ele=2; % order of filter for EEG data
Settings.HP_acc = 1; % Hz
Settings.order_acc=2; % order of filter for accelerometer data
Settings.normFreq_ele=50; % Hz upper limit for normalized spectrogram
Settings.Plot=true;
Settings.SelectMultipleChannels=true;
Settings.ShowLinearFit=false;
Settings.Randomize=false;
Settings.time_threshold=0;

Channels={'right-S';'right-M';'left-M';'left-S'};
nChannels=length(Channels);

% Colors for plotting
c=[0,0.4470,0.7410;0.8500,0.3250,0.0980;0.4660,0.6740,0.1880]; % c=[blue;red;green]

%% Select to be analyzed .mat files
[E_Name,PathName] = uigetfile('*.mat','Select the file to analyse','MultiSelect', 'on');
if ischar(E_Name) %if-statement for the case that only one file is selected
    E_Name={E_Name}; 
end
cd(PathName)
nE=length(E_Name);

% Assume that the sampling rate is constant across measurements
load(E_Name{1},'Electrical')
Settings.fs=Electrical.fs;

%% Import experiment data
% import data from .mat file
select_points=cell(nE,1);
E_cellTable=cell(nE,1);
PSD=cell(nE,nChannels+1);
subplotIndex=reshape(1:15,[3,5])';
for i=1:nE
    
    load(E_Name{i},'Acceleration','Electrical')
    Acceleration.t=(0:1/Acceleration.fs:(size(Acceleration.XYZ,1)-1)/Acceleration.fs)'; % s
    Electrical.t=(0:1/Electrical.fs:(size(Electrical.CH1234,1)-1)/Electrical.fs)'; % s
    
    %% Approximate dynamic acceleration by applying a high-pass filter and calculate length of vector
    % approximate XYZ dynamic acceleration by applying a high-pass filter
    Acceleration.dyn=[HighPassfilter(Settings.order_acc,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,1)),...
                      HighPassfilter(Settings.order_acc,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,2)),...
                      HighPassfilter(Settings.order_acc,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,3))];
    Acceleration=rmfield(Acceleration,'XYZ');
                  
    % calculate length of XYZ dynamic acceleration vector
    Acceleration.dyn=vecnorm(Acceleration.dyn,2,2);

    % interpolate acceleration to same time grid as electrical
    Acceleration.dyn=interpn(Acceleration.t,Acceleration.dyn,Electrical.t);
    
    % select part of trace
    Settings.fs=Settings.fs/5;
    select_points{i} = EEG_Select_E([Electrical.CH1234(1:5:end,:),Acceleration.dyn(1:5:end)],Electrical.t(1:5:end),[],E_Name{i},Settings,[Channels;'Accelerometer']);
    close(gcf)
    for ii=1:nChannels+1
        select_points{i}{ii} = 5*select_points{i}{ii}; 
    end
    Settings.fs=Settings.fs*5;
    
    % Filter
    for ii=1:nChannels
        Electrical.CH1234(:,ii) = HighLowPassfilter(Settings.order_ele,...
        [Settings.HP_ele,Settings.LP_ele], Electrical.fs, Electrical.CH1234(:,ii));
    end
        
    % Initialize figure, plots added under "compute PSD"
    figure('name',E_Name{i})
    
    % compute PSD
    nWindowLength=round(Settings.window*Settings.fs);
    nOverlap=round(nWindowLength*Settings.overlap);
    nFFT=nWindowLength;
    for ii=1:nChannels+1
        if ii<nChannels+1
            [PSD_matrix,t,f] = spectrogram_dis(Electrical.CH1234(:,ii),select_points{i}{ii},nWindowLength,nOverlap,nFFT,Settings.fs);
            PSD{i,ii}=mean(PSD_matrix,2);
        else
            [PSD_matrix,t,f] = spectrogram_dis(Acceleration.dyn,select_points{i}{ii},nWindowLength,nOverlap,nFFT,Settings.fs);
        end
        
        if Settings.Plot
            % Trace
            subplot(nChannels+1,3,subplotIndex(ii,1)), hold on
            temp_t=Electrical.t(1:5:end); 
            if ii<nChannels+1
                temp_y=Electrical.CH1234(1:5:end,ii);
                plot(temp_t,temp_y)
                select_points_temp=select_points{i}{ii}/5;
                for iii=1:size(select_points{i}{ii},1)
                    plot(temp_t(select_points_temp(iii,1):select_points_temp(iii,2)),...
                         temp_y(select_points_temp(iii,1):select_points_temp(iii,2)),...
                         'color',c(3,:))
                end
                xlabel('Time (s)')
                title(Channels{ii},'Interpreter','none')
                ylabel('EEG (\muV)')
            else
                temp_y=Acceleration.dyn(1:5:end);
                plot(temp_t,temp_y)
                select_points_temp=select_points{i}{ii}/5;
                for iii=1:size(select_points{i}{ii},1)
                    plot(temp_t(select_points_temp(iii,1):select_points_temp(iii,2)),...
                         temp_y(select_points_temp(iii,1):select_points_temp(iii,2)),...
                         'color',c(3,:))
                end
                xlabel('Time (s)')
                title('Accelerometer','Interpreter','none')
                ylabel('Dynamic acceleration (g)')
            end
            xlim(Electrical.t([1,end])')
            % Spectrogram    
            subplot(nChannels+1,3,subplotIndex(ii,2))
            t_spec=[0,sum(select_points{i}{ii}(:,2)-select_points{i}{ii}(:,1)+1)/Settings.fs];
            imagesc(t_spec,f,PSD_matrix)
            if ii<nChannels+1
                title(Channels{ii},'Interpreter','none')
            else
                title('Accelerometer','Interpreter','none')
            end
            xlabel('Time (s)'),ylabel('Frequency (Hz)')
            set(gca,'YDir','normal')
            set(gca,'CLim',[0,mean(max(PSD_matrix))])
            ylim([0,20])
            % PSDs
            if ii<nChannels+1
                subplot(nChannels+1,3,subplotIndex(:,3))
                plot(f,PSD{i,ii}),hold on
            end
        end
    end
    subplot(nChannels+1,3,subplotIndex(:,3))
    xlabel('Frequency (Hz)'),ylabel('PSD (\muV^2)')
    title('Spectra of each EEG channel','Interpreter','none')
    xlim([0,Settings.normFreq_ele])
    legend(Channels)
    grid on
    
    % Create cellTable for each experiment
    E_cellTable{i}=cell(length(f)+2,nChannels+1);
    E_cellTable{i}(1,1)=E_Name(i);
    for ii=1:nChannels
        E_cellTable{i}(2,ii)={[Channels{ii},' (muV^2)']};
        E_cellTable{i}(3:end,ii)=num2cell(PSD{i,ii});
    end
    
end

%% Save2Excel
charDateTime=char(datetime);
charDateTime(charDateTime==':')='-';
xlswrite(['SpectraSelectionAnalysis ',charDateTime,'.xlsx'],[[[{[]};{'Frequency (Hz)'};num2cell(f)],cell(length(f)+2,1)],cat(2,E_cellTable{:})])