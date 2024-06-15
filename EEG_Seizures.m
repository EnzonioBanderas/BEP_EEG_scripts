% Plot PSD and spectrogram of acceleration and electrical data for each
% selected file
clearvars, close all
%% Input
% Settings
Settings.t_prepost=200; % s
Settings.window=2; % s
Settings.overlap_ratio=0.5; % ratio
Settings.HP_ele = 1; % Hz
Settings.LP_ele = 100; % Hz
Settings.HP_acc = 1; % Hz
Settings.maxFreq_acc=50; % Hz upper limit for normalized spectrogram
Settings.maxFreq_ele=50; % Hz upper limit for normalized spectrogram

% Select to be analyzed .mat files
[FileName,PathName] = uigetfile('*.mat','Select the file to analyse','MultiSelect', 'on');
if ischar(FileName) %if-statement for the case that only one file is selected
    FileName={FileName}; 
end
cd(PathName)
nExperiments=length(FileName);

%% Plot setttings
% plot property 'Color'
Settings.c=colormap(lines); % colors to most contrasting color map
close all

% plot property 'LineStyle'
Settings.ls=["-";"--";":";"-."];

%% connect file names with seizure times
% "Manual" Input, can be altered as (standard) input window
% seizure_startend=zeros(nExperiments,2);
% for i=1:nExperiments
%     switch FileName{i}(1:end-4)
%         case '10425-01_het_D1'
%             seizure_startend(i,:)=[170417.285,170590.522];
%         case '10425-01_het_D1'
%             seizure_startend(i,:)=[170417.285,170590.522];
%         case '10425-01_het_D1'
%             seizure_startend(i,:)=[170417.285,170590.522];
%         case '10425-01_het_D1'
%             seizure_startend(i,:)=[170417.285,170590.522];
%         case '10425-01_het_D1'
%             seizure_startend(i,:)=[170417.285,170590.522];
%         case '10425-01_het_D1'
%             seizure_startend(i,:)=[170417.285,170590.522];
%         case '10425-01_het_D1'
%             seizure_startend(i,:)=[170417.285,170590.522];
%         case '10425-01_het_D1'
%             seizure_startend(i,:)=[170417.285,170590.522];
%     end
% end

% Old loop
seizure_start=zeros(nExperiments,1);
for i=1:nExperiments
    switch FileName{i}(1:end-4)
        case '10425-01_het_D1'
            seizure_start(i)=171006.069;
        case '10425-02_wt_D1'
            seizure_start(i)=170681.543;
        case '10425-04_het_D1'
            seizure_start(i)=170858.569;
        case '10425-05_wt_D1'
            seizure_start(i)=170973.901;
        case '10425-05_wt_D2'
            seizure_start(i)=175012.231;
        case '10425-06_het_D1'
            seizure_start(i)=NaN; % no signal due to loss of reference
        case '10430-02_het_D1'
            seizure_start(i)=175094.019;
        case '10430-05_het_D1'
            seizure_start(i)=175277.388;
        case '10630-01_wt_D1'
            seizure_start(i)=175359.391;
        case '10630-03_wt_D1'
            seizure_start(i)=175646.235;
        case '10631-02_wt_D1'
            seizure_start(i)=NaN; % bad trace, will be repeated
    end
end

%% Start of experiment loop
for i=1:nExperiments
    
%% Import data
load(FileName{i},'Acceleration','Electrical')

%% Pre processing (time vectors, seizure logical, truncation)
% create time vectors
Acceleration.t=0:1/Acceleration.fs:(size(Acceleration.XYZ,1)-1)/Acceleration.fs; % s
Acceleration.nDim=size(Acceleration.XYZ,2);
Electrical.t=0:1/Electrical.fs:(size(Electrical.CH1234,1)-1)/Electrical.fs; % s
Electrical.nChannel=size(Electrical.CH1234,2);

% create logical
Acceleration.t_logical=Acceleration.t>seizure_start(i)-Settings.t_prepost&...
                       Acceleration.t<seizure_start(i)+Settings.t_prepost;
Electrical.t_logical=Electrical.t>seizure_start(i)-Settings.t_prepost&...
                     Electrical.t<seizure_start(i)+Settings.t_prepost;
                   
% truncation
Acceleration.XYZ=Acceleration.XYZ(Acceleration.t_logical,:);
Acceleration.t=Acceleration.t(Acceleration.t_logical);
Electrical.CH1234=Electrical.CH1234(Electrical.t_logical,:);
Electrical.t=Electrical.t(Electrical.t_logical);

%% Plot and calculation, first 3 Electrical figures, second 4 Acceleration figures
% Electrical
% 1) (filtered) normalized spectrogram
% 2) PSD
% 3) pre and post PSDs
% Acceleration
% 1) dynamic acceleration vector components
% 2) dynamic acceleration vector length
% 3) unfiltered normalized spectrogram
% 4) filtered normalized spectrogram

% FFT pre-processing
temp.window = round( Electrical.fs * Settings.window );
temp.noverlap = round( temp.window * Settings.overlap_ratio );
temp.nfft = temp.window;

    % 3 Figures
    fig.trace=figure('Name','Filtered electrical trace');
    fig.normspec=figure('Name','Normalized spectrograms');
    fig.PSD=figure('Name','PSD');
    fig.prepost_PSD=figure('Name','prepost_PSD');
    for ii=1:Electrical.nChannel % Loop through channels
    
    % Filter Electrical signals      
    Electrical.CH1234(:,ii) = HighPassfilter(2, Settings.HP_ele, Electrical.fs, Electrical.CH1234(:,ii));
    Electrical.CH1234(:,ii) = LowPassfilter(2, Settings.LP_ele, Electrical.fs, Electrical.CH1234(:,ii));
    
    % Plot Filtered Electrical signals
    set(0,'CurrentFigure',fig.trace)
    subplot(Electrical.nChannel,1,ii)
    temp.ttt=Electrical.t-Electrical.t(1)-Settings.t_prepost;
    plot(temp.ttt,Electrical.CH1234(:,ii))
    xlabel('Time (s)'),ylabel('Amplitude (\muV)')
    title(['CH',num2str(ii)])
    
    % STFT
    [~,temp.fff,temp.ttt,temp.ps]=spectrogram(Electrical.CH1234(:,ii),...
    temp.window, temp.noverlap, temp.nfft, Electrical.fs);

    % Adjust time
    temp.ttt=temp.ttt-Settings.t_prepost;

    % Normalization
    temp.ps_norm=(temp.ps./bandpower2(temp.ps,temp.fff,[Settings.HP_ele,Settings.maxFreq_ele]))*100; %normalized power spectrum with frequencies of interest (%)

    % Normalized (Settings.HP_ele-maxFreq Hz) spectrogram plot
    set(0,'CurrentFigure',fig.normspec)
    subplot(Electrical.nChannel,1,ii)
    imagesc(temp.ttt,temp.fff,temp.ps_norm)
    set(gca,'YDir','normal')
    colormap('PARULA')
    caxis([0,mean(max(temp.ps_norm))])
    h=colorbar;
    h.Label.String = 'PSD (%)';  
    xlabel('Time (s)'),ylabel('Frequency (HZ)')
    title(['CH',num2str(ii)])

    % PSD calculation
    temp.PSD=mean(temp.ps,2);

    % PSD plot
    set(0,'CurrentFigure',fig.PSD)
    subplot(Electrical.nChannel,1,ii)
    plot(temp.fff,temp.PSD)
    xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    title(['CH',num2str(ii)])

    % pre-PSD calculation
    temp.logical_pre=temp.ttt<0;
    temp.PSD_pre=mean(temp.ps(:,temp.logical_pre),2);
    % post-PSD calculation
    temp.logical_post=temp.ttt>0;
    temp.PSD_post=mean(temp.ps(:,temp.logical_post),2);
    % prepost PSD plot
    set(0,'CurrentFigure',fig.prepost_PSD)
    % pre
    subplot(Electrical.nChannel,2,(ii-1)*2+1)
    plot(temp.fff,temp.PSD_pre)
    xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    title(['CH',num2str(ii),' -pre'])
    % post
    subplot(Electrical.nChannel,2,ii*2)
    plot(temp.fff,temp.PSD_post)
    xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    title(['CH',num2str(ii),' -post'])

    end
    
    % approximate XYZ dynamic acceleration by applying a high-pass filter
    Acceleration.dyn_XYZ=[HighPassfilter(2,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,1)),...
                          HighPassfilter(2,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,2)),...
                          HighPassfilter(2,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,3))];
                  
    % calculate length of XYZ dynamic acceleration vector
    Acceleration.dyn_length=vecnorm(Acceleration.dyn_XYZ,2,2);
                  
    % dynamic acceleration vector length plot
    fig.dynLength=figure('Name','dynamical acceleration vector length over time around seizure');
    plot(Acceleration.t-seizure_start(i),Acceleration.dyn_length);
    xlabel('Time (s)'), ylabel('Dynamical acceleration vector length (g -probably)')
    
    % dynamic acceleration vector components plot
    fig.dynXYZ=figure('Name','dynamical acceleration vector XYZ components over time around seizure');
    xlabel('Time (s)'), ylabel('Dynamical acceleration (g)')
    
    % FFT pre-processing
    temp.window = round( Acceleration.fs * Settings.window );
    temp.noverlap = round( temp.window * Settings.overlap_ratio );
    temp.nfft = temp.window;
    
    % normalized spectrogram plots unfiltered and filtered
    fig.dynSpec_unfilt=figure('Name','dynamical acceleration spectrogram of time around seizure -unfiltered');
    fig.dynSpec=figure('Name','dynamical acceleration spectrogram of time around seizure');
    
    for ii=1:Acceleration.nDim
                   
    % plot XYZ components
    set(0,'CurrentFigure',fig.dynXYZ)
    plot(Acceleration.t-seizure_start(i),Acceleration.dyn_XYZ(:,ii),'Color',Settings.c(ii,:)), hold on
    
    %%%%%%%%%%
    % spec calculations
    [~,temp.fff,temp.ttt,temp.ps]=spectrogram(Acceleration.XYZ(:,ii),...
    temp.window, temp.noverlap, temp.nfft, Acceleration.fs);

    % Adjust time
    temp.ttt=temp.ttt-Settings.t_prepost;

    % plot normalized (0-maxFreq Hz) spectrogram
    set(0,'CurrentFigure',fig.dynSpec_unfilt)
    subplot(Acceleration.nDim,1,ii)  
    temp.logical_fff=temp.fff<=Settings.maxFreq_acc;
    temp.ps_norm=(temp.ps(temp.logical_fff,:)./sum(temp.ps(temp.logical_fff,:)))*100; %normalized power spectrum with frequencies of interest (%)
    imagesc(temp.ttt,temp.fff(temp.logical_fff),temp.ps_norm)
    set(gca,'YDir','normal')
    colormap('PARULA')
    caxis([0,mean(max(temp.ps_norm))])
    h=colorbar;
    h.Label.String = 'PSD (%)';
    xlabel('Time (s)'),ylabel('Frequency (HZ)')
%     ,clabel('Acceleration PSD (%)')
    
    switch ii
        case 1
        title("dynamical acceleration X component");
        case 2
        title("dynamical acceleration Y component");
        case 3
        title("dynamical acceleration Z component");
    end
    %%%%%%%
    
    %%%%%%%%%%
    % spec calculations
    [~,temp.fff,temp.ttt,temp.ps]=spectrogram(Acceleration.dyn_XYZ(:,ii),...
    temp.window, temp.noverlap, temp.nfft, Acceleration.fs);

    % Adjust time
    temp.ttt=temp.ttt-Settings.t_prepost;

    % plot normalized (0-maxFreq Hz) spectrogram
    set(0,'CurrentFigure',fig.dynSpec)
    subplot(Acceleration.nDim,1,ii)  
    temp.logical_fff=temp.fff<=Settings.maxFreq_acc&temp.fff>=Settings.HP_acc;
    temp.ps_norm=(temp.ps(temp.logical_fff,:)./sum(temp.ps(temp.logical_fff,:)))*100; %normalized power spectrum with frequencies of interest (%)
    imagesc(temp.ttt,temp.fff(temp.logical_fff),temp.ps_norm)
    set(gca,'YDir','normal')
    colormap('PARULA')
    caxis([0,mean(max(temp.ps_norm))])
    h=colorbar;
    h.Label.String = 'PSD (%)';
    xlabel('Time (s)'),ylabel('Frequency (HZ)')
%     ,clabel('Acceleration PSD (%)')
    
    switch ii
        case 1
        title("dynamical acceleration X component");
        case 2
        title("dynamical acceleration Y component");
        case 3
        title("dynamical acceleration Z component");
    end
    %%%%%%%%%%%%%%
    
    end % end of acceleration dimension loop

% clearvars -except FileName nExperiments
end

