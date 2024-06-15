clearvars, close all
%% Input 
Settings.window=2; % s
Settings.overlap_ratio=0.5; % ratio
Settings.HP_ele = 1; % Hz
Settings.LP_ele = 100; % Hz
Settings.HP_acc = 1; % Hz
Settings.p_value=.03; % probability
Settings.normFreq_ele=50; % Hz upper limit for normalized spectrogram
Settings.timeBin=Settings.window*30; % s
Settings.overlapBin=0.5; % ratio
Settings.timeSelect=4*60*60; % s
Settings.StateString=["Awake","Sleep"]; % states that are to be selected for (Awake/Sleep)

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
Electrical.t=(0:1/Electrical.fs:(size(Electrical.CH1234,1)-1)/Electrical.fs)'; % s
Acceleration.t=(0:1/Acceleration.fs:(size(Acceleration.XYZ,1)-1)/Acceleration.fs)'; % s

% truncate acceleration % can still be changed to be correct for overlap
% with buffer overlap truncation becomes unnecessary
Acceleration.Select=DataTable_Logical.Select_Acceleration{temp.select_index};
Acceleration.Select_points=logical2points(Acceleration.Select);
temp.segment_length=Acceleration.Select_points(:,2)-Acceleration.Select_points(:,1)+1;
Acceleration.nBin=round(Settings.timeBin*Acceleration.fs);
Acceleration.Select_points=Acceleration.Select_points(:,2)-mod(temp.segment_length,Acceleration.nBin);
Acceleration.Select=points2logical(Acceleration.Select_points);
Acceleration.XYZ=Acceleration.XYZ(Acceleration.Select,:);
Acceleration.t=Acceleration.t(Acceleration.Select);

% truncate electrical
Electrical.Select=DataTable_Logical.Select_E;ectrical{temp.select_index};
Electrical.Select_points=logical2points(Electrical.Select);
temp.segment_length=Electrical.Select_points(:,2)-Electrical.Select_points(:,1)+1;
Electrical.nBin=round(Settings.timeBin*Electrical.fs);
Electrical.Select_points=Electrical.Select_points(:,2)-mod(temp.segment_length,Electrical.nBin);
Electrical.Select=points2logical(Electrical.Select_points);
Electrical.CH1234=Electrical.CH1234(Electrical.Select,:);
Electrical.t=Electrical.t(Electrical.Select);

% some constants for acceleration data
Acceleration.n=length(Acceleration.t);
Acceleration.nCol=Acceleration.n/Acceleration.nBin;
Acceleration.nDim=size(Acceleration.XYZ,2);
% some constants for electrical data
Electrical.n=length(Electrical.t);
Electrical.nCol=Electrical.n/Electrical.nBin;
Electrical.nChannel=size(Electrical.CH1234,2);

clearvars temp
%% Approximate dynamic acceleration by applying a high-pass filter and calculate length of vector
% approximate XYZ dynamic acceleration by applying a high-pass filter
Acceleration.dyn=[HighPassfilter(2,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,1)),...
                  HighPassfilter(2,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,2)),...
                  HighPassfilter(2,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,3))];
                   
% calculate length of XYZ dynamic acceleration vector
Acceleration.dyn=vecnorm(Acceleration.dyn,2,2);

%% Select 
% Get points of consecutive sequences of time trace for selection which has
% a length defined by Settings.timeSelect (s)
temp.nRow_acc=round(Settings.timeSelect*Acceleration.fs);
temp.nRow_ele=round(Settings.timeSelect*Electrical.fs);
temp.nCol_acc=floor(Acceleration.n/temp.nRow_acc);
temp.nCol_ele=floor(Electrical.n/temp.nRow_ele);
[~,~,temp.points_acc]=reshape2(Acceleration.t,[temp.nRow_acc,temp.nCol_acc]);
[~,~,temp.points_ele]=reshape2(Electrical.t,[temp.nRow_ele,temp.nCol_ele]);
       
% Create figure       
figure('Name',[Mouse,' Experiment number ',num2str(experiment_index)],...
       'units','normalized','outerposition',[0 0 1 1])

% preallocate false values to state logical for acceleraion end electrical data
select_logical_acc=cell(nState,1);
for i=1:nState
    select_logical_acc{i}=false(Acceleration.n,1);
end
select_logical_ele=cell(nState,1);
for i=1:nState
    select_logical_ele{i}=false(Electrical.n,1);
end

% select points and assign logicals for each state for acceleration and electrical
for j=1:size(temp.points_acc,1)
    
    time_vector_unselec_acc=(Acceleration.t(temp.points_acc(j,1):temp.points_acc(j,2)))';
    time_vector_unselec_ele=(Electrical.t(temp.points_ele(j,1):temp.points_ele(j,2)))';
    
    dyn_vector_unselec=Acceleration.dyn(temp.points_acc(j,1):temp.points_acc(j,2));
    
    plot(time_vector_unselec_acc,dyn_vector_unselec,'Color',Settings.c(nState+1,:))
    hold on
    xlabel('Time (s)'),ylabel('Acceleration (g)')
    xlim(time_vector_unselec_acc([1,end]))
    title(['Part ',num2str(j),' ',num2str(Acceleration.t(temp.points_acc(j,1))),'-',num2str(Acceleration.t(temp.points_acc(j,2))),' s'])
    grid on
    
        % Select time frame
        for jj=1:nState
        
        m=msgbox(['Select ',Settings.StateString{jj},' time frames?']);
        waitfor(m)
        
        pointx_temp=0; index=0;
        while ~isempty(pointx_temp)
            index=index+1;
            [pointx_temp,~]=ginput(1);
            if ~isempty(pointx_temp)
                
                pointx(index)=pointx_temp;
                plot(pointx(index),dyn_vector_unselec(dsearchn(time_vector_unselec_acc,pointx(index))),'rx')
            
                if mod(index,2)==0
                    % assign values to temporary acceleration logical
                    select_logical_temp=time_vector_unselec_acc>=pointx(index-1)&...
                                        time_vector_unselec_acc<=pointx(index);
                    select_logical_acc{jj}(temp.points_acc(j,1):temp.points_acc(j,2))=...
                    select_logical_acc{jj}(temp.points_acc(j,1):temp.points_acc(j,2))|select_logical_temp;
                    
                    % plot
                    plot(time_vector_unselec_acc(select_logical_temp),...
                         dyn_vector_unselec(select_logical_temp),...
                         'Color',Settings.c(jj,:))
                    
                    % assign values to temporary electrical logical
                    select_logical_temp=time_vector_unselec_ele>=pointx(index-1)&...
                                        time_vector_unselec_ele<=pointx(index);
                    select_logical_ele{jj}(temp.points_ele(j,1):temp.points_ele(j,2))=...
                    select_logical_ele{jj}(temp.points_ele(j,1):temp.points_ele(j,2))|select_logical_temp;
                end
            
            end
        end
        
        end
        cla % clear axes
    
end
close(gcf)

% Assign state logicals to Acceleration and Electrical structures
for jj=1:nState
eval(['Acceleration.',Settings.StateString{jj},'=select_logical_acc{jj};'])
eval(['Electrical.',Settings.StateString{jj},'=select_logical_ele{jj};'])
end

Acceleration.indeterminate=sum([select_logical_acc{:}],2)==0;
Electrical.indeterminate=sum([select_logical_ele{:}],2)==0;

% Assign indeterminate
Acceleration.indeterminate=~(Acceleration.Awake|Acceleration.Sleep);
Electrical.indeterminate=~(Electrical.Awake|Electrical.Sleep);

% Convert logical arrays to points
Acceleration.Awake_points=logical2points(Acceleration.Awake);
Electrical.Awake_points=logical2points(Electrical.Awake);
Acceleration.Sleep_points=logical2points(Acceleration.Sleep);
Electrical.Sleep_points=logical2points(Electrical.Sleep);
Acceleration.indeterminate_points=logical2points(Acceleration.indeterminate);
Electrical.indeterminate_points=logical2points(Electrical.indeterminate);

clearvars temp
%% Filter Electrical signals        
for j=1:Electrical.nChannel % Loop through channels
Electrical.CH1234(:,j) = HighPassfilter(2, Settings.HP_ele, Electrical.fs, Electrical.CH1234(:,j));
Electrical.CH1234(:,j) = LowPassfilter(2, Settings.LP_ele, Electrical.fs, Electrical.CH1234(:,j));
end

%% Start of plotting
%% Plot of electrical data
temp.points={Electrical.Awake_points,Electrical.Sleep_points,Electrical.indeterminate_points};
figure('Name','')
for j=1:Electrical.nChannel
    subplot(Electrical.nChannel,1,j)
    hold on
    xlabel('Time (s)'),ylabel('Amplitude (mV)')
    title(['CH',num2str(j)])
    xlim(Electrical.t([1,end]))
    for jj=1:length(temp.points)
        for jjj=1:size(temp.points{jj},1)
            temp.index=temp.points{jj}(jjj,1):temp.points{jj}(jjj,2);
            plot(Electrical.t(temp.index),Electrical.CH1234(temp.index,j),...
                 'Color',Settings.c(jj,:))
        end
    end
end
%% Length of dynamical acceleration vector plot
temp.points={Acceleration.Awake_points,Acceleration.Sleep_points,Acceleration.indeterminate_points};

figure('Name',[Mouse,' Change of acceleration state plot'])
for j=1:nState+1 % state
    for jj=1:size(temp.points{j},1)
        
    plot(Acceleration.t(temp.points{j}(jj,1):temp.points{j}(jj,2)),...
         Acceleration.dyn(temp.points{j}(jj,1):temp.points{j}(jj,2)),...
         'Color',Settings.c(j,:)), hold on

    end
end
xlabel('Time (s)'), ylabel('Dynamic acceleration (g)')

clearvars temp
%% Calculate and plot FFT's
% FFT pre-processing
temp.window =   round(Electrical.fs*Settings.window);
temp.noverlap = round(temp.window*Settings.overlap_ratio);
temp.nfft = temp.window;

temp.points={Electrical.Awake_points,Electrical.Sleep_points};
Electrical.PSD_cell=cell(2,Electrical.nChannel);

% calculate PSD with spectrogram with points
figure('Name',[Mouse,' Normalized Spectrogram'])
for j=1:nState % State
for jj=1:Electrical.nChannel % Channel
    
    temp.nPoints=size(temp.points{j},1);
    temp.ps=cell(temp.nPoints,1);
    for jjj=1:temp.nPoints

    [~,~,~,temp.ps{jjj}]=spectrogram(Electrical.CH1234(temp.points{j}(jjj,1):temp.points{j}(jjj,2),jj),...
    temp.window, temp.noverlap, temp.nfft, Electrical.fs);

    end
    temp.ps=cat(2,temp.ps{:});
    
% plot normalized (Settings.HP_ele-maxFreq Hz) spectrogram
subplot(Electrical.nChannel,2,j+(jj-1)*2)
temp.fff=0:Electrical.fs/temp.nfft:Electrical.fs/2; % approximate frequency bins positions directly
temp.ttt=Settings.window*Settings.overlap_ratio:Settings.window*Settings.overlap_ratio:size(temp.ps,2); % approximate center time position of PSDs directly
temp.ps_norm=(temp.ps./bandpower2(temp.ps,temp.fff,[Settings.HP_ele,Settings.normFreq_ele]))*100; %normalized power spectrum with frequencies of interest (%)
imagesc(temp.ttt,temp.fff,temp.ps_norm)
set(gca,'YDir','normal')
colormap('PARULA')
caxis([0,mean(max(temp.ps_norm))])
colorbar              
xlabel('Time (s)'),ylabel('Lfs PSD (%)')
ylim([Settings.HP_ele,Settings.normFreq_ele])

if j==1
    title(['Awake CH',num2str(jj)])
else
    title(['Sleep CH',num2str(jj)]);
end

Electrical.PSD_cell{j,jj}=mean(temp.ps,2); % rows for Awake/Sleep, columns for channels

end
end

%% PSD plot
figure('Name',[Mouse,' PSD'])
for j=1:nState % state
for jj=1:Electrical.nChannel % channel
    
subplot(Electrical.nChannel,2,j+(jj-1)*2)
plot(temp.fff,Electrical.PSD_cell{j,jj})
xlim([Settings.HP_ele,50])
    
if j==1
    title(['Awake CH',num2str(jj)])
else
    title(['Sleep CH',num2str(jj)]);
end
    
end
end
xlabel('Frequency (Hz)'), ylabel('Power \muV^2')

clearvars temp
%% Convert variables to DataTable 
Mouse=strsplit(Mouse,'_');
Name=Mouse(1);
Genotype=Mouse(2); 
Day=Mouse(3);
Channel=["CH1";"CH2";"CH3";"CH4"];

nPSD_cell=numel(Electrical.PSD_cell);

Name_string=strings(nPSD_cell,1);
Genotype_string=Name_string;
Day_string=Name_string;
Channel_string=Name_string;
State_string=Name_string;
PSD_table_cell=cell(nPSD_cell,1);
Electrical_fs=zeros(nPSD_cell,1);

Name_string(1:end)=Name;
Genotype_string(1:end)=Genotype;
Day_string(1:end)=Day;

count=0;
for d=1:nState % States
    for dd=1:Electrical.nChannel % Channels
    
    count=count+1;    
        
    Channel_string(count)=Channel(dd);
    State_string(count)=Settings.StateString(d);
    PSD_table_cell(count)=Electrical.PSD_cell(d,dd);
    Electrical_fs(count)=Electrical.fs;
    eval(['Acceleration_Points_cell{count}=Acceleration.',State_string(count),'_points;'])
    eval(['Electrical_Points_cell{count}=Electrical.',State_string(count),'_points;'])
    Truncation_Points_cell{count}=Truncation_Points;

    end
end

TableNames={'Name','Genotype','Day','Channel','State','Power','fs','Truncation_Points','Acceleration_Points','Sleep_Points'};
DataTable_cell{experiment_index}=table(Name_string,Genotype_string,Day_string,Channel_string,State_string,...
                PSD_table_cell,Electrical_fs,Truncation_Points_cell,Acceleration_Points_cell,Electrical_Points_cell,...
                'VariableNames',TableNames);

clearvars -except Settings FileName PathName nExperiments experiment_index DataTable_cell DataTable_Logical nState
end

% save datatable and other settings
cd ..
DataTable=cat(1,DataTable_cell{:});
save(['DataTable_Manual_',num2str(now),'.mat'],'DataTable')
save('Settings','Settings')