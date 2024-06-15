%% Deadline 02-11-2018
%
% For a continuous segment that is >=1h, split it into 0.5h overlapping
% segments of 1h. Simply include segments that are >0.5h&<1h without
% splitting. Calculate percentage in sleep and percentage REM in sleep 
% (with this you can get all other percentages). 
%
% Plot barplots AS in correct order.
%
% Add "sigstar" to AS barplots programmatically.

% compile *_Analyzed.mat files in folder into DataTable's which contain PSD
% and state information
clearvars,close all
%% Get folder
EF = uigetdir('','Select Experiment Folder (EF)');
cd(EF)
E=dir('*_Analyzed.mat');
nE=length(E);

% Load start time table
load('DataTable_StartTime')
load('Settings')

% Define between which hours it is light and dark
Light=[8,18];
Dark=[20,6];

% preallocation
DataTable2_cell=cell(nE,1);
DataTable3_cell=cell(nE,1);
DataTable4_cell=cell(nE,1);
DataTable5_cell=cell(nE,1);
DataTable_Percentages_cell=cell(nE,1);
Name_Percentages_cell=cell(nE,1);
Genotype_Percentages_cell=cell(nE,1);
Day_Percentages_cell=cell(nE,1);
Light_Percentages_cell=cell(nE,1);
Dark_Percentages_cell=cell(nE,1);

% DataTable6 preallocation
logseg_states={'Awake';'NREM';'REM'};
logseg_nState=length(logseg_states);
D6_Name    =cell(nE*logseg_nState,1);
D6_Genotype=cell(nE*logseg_nState,1);
D6_Day     =cell(nE*logseg_nState,1);
D6_State   =cell(nE*logseg_nState,1);
D6_Data    =cell(nE*logseg_nState,1);

%% Light/Dark assume to be from 7 AM to 7 PM, compute percentages and state length distributions
bar = waitbar(0, 'Compiling Data, please wait');
for i=1:nE
waitbar((i-1)/nE, bar, 'Compiling Data, please wait');

% States
States={'Awake_Light';'NREM_Light';'REM_Light';...
        'Awake_Dark';'NREM_Dark';'REM_Dark'};
nState=length(States);

% Suggestion for selection
cd('C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Clean')
load([E(i).name(1:end-13),'.mat'],'Electrical')
suggestion=true(size(Electrical.CH1234,1),1);
kernel=true(round(20*Electrical.fs),1); % Kernel of 5 s
for ii=1:Settings.nChannel
    suggestion=suggestion&~conv(abs(Electrical.CH1234(:,ii))>max(Electrical.CH1234(:,ii))-100,kernel,'same');
end
cd(EF)

    % Load and create Light and Dark logicals
    load(E(i).name,'Electrical','FFT','Settings')
    startTime=DataTable_StartTime.StartTime(strcmp(DataTable_StartTime.Name,E(i).name));
    startTime=hour(startTime)*60*60+minute(startTime)*60+second(startTime);
    t=mod(Electrical.t+startTime,24*60*60);
    Light_logical=t>Light(1)*60*60&t<Light(2)*60*60;
    Dark_logical=t>Dark(1)*60*60|t<Dark(2)*60*60;
    Light_logical=Light_logical&points2logical(Electrical.points,Electrical.n);
    Dark_logical=Dark_logical&points2logical(Electrical.points,Electrical.n);
    points_Light=logical2points(Light_logical);
    points_Dark=logical2points(Dark_logical);
    
    % Temptable strings
    Esplit=strsplit(E(i).name,'_');
    E_Name=Esplit{1};
    E_Genotype=Esplit{2};
    E_Day=Esplit{3};

    % Distributions
    Awake_logical=points2logical(Electrical.points_AA,Electrical.n)|points2logical(Electrical.points_QA,Electrical.n);
    NREM_logical =points2logical(Electrical.points_NREM,Electrical.n);
    REM_logical  =points2logical(Electrical.points_REM,Electrical.n);
    % Light
    Electrical_points.Awake_Light=logical2points(Awake_logical&Light_logical);
    Awake_Light=(Electrical_points.Awake_Light(:,2)-Electrical_points.Awake_Light(:,1)+1)/Electrical.fs;
    Electrical_points.NREM_Light=logical2points(NREM_logical&Light_logical);
    NREM_Light=(Electrical_points.NREM_Light(:,2)-Electrical_points.NREM_Light(:,1)+1)/Electrical.fs;
    Electrical_points.REM_Light=logical2points(REM_logical&Light_logical);
    REM_Light=(Electrical_points.REM_Light(:,2)-Electrical_points.REM_Light(:,1)+1)/Electrical.fs;
    % Dark
    Electrical_points.Awake_Dark=logical2points(Awake_logical&Dark_logical);
    Awake_Dark=(Electrical_points.Awake_Dark(:,2)-Electrical_points.Awake_Dark(:,1)+1)/Electrical.fs;
    Electrical_points.NREM_Dark=logical2points(NREM_logical&Dark_logical);
    NREM_Dark=(Electrical_points.NREM_Dark(:,2)-Electrical_points.NREM_Dark(:,1)+1)/Electrical.fs;
    Electrical_points.REM_Dark=logical2points(REM_logical&Dark_logical);
    REM_Dark=(Electrical_points.REM_Dark(:,2)-Electrical_points.REM_Dark(:,1)+1)/Electrical.fs;
    % Percentages Light and Dark
    Light_Percentages=[sum(Awake_Light),sum(NREM_Light),sum(REM_Light)];
    Light_Percentages(isempty(Light_Percentages))=0;
    Light_Percentages=Light_Percentages/sum(Light_Percentages);
    Light_Percentages(3)=Light_Percentages(3)/sum(Light_Percentages([2,3]));
    Dark_Percentages=[sum(Awake_Dark),sum(NREM_Dark),sum(REM_Dark)];
    Dark_Percentages(isempty(Dark_Percentages))=0;
    Dark_Percentages=Dark_Percentages/sum(Dark_Percentages);
    Dark_Percentages(3)=Dark_Percentages(3)/sum(Dark_Percentages([2,3]));
    
    % DataTable6
    segmentLength=round(60*60*Electrical.fs);
    [Awake_logseg,Awake_midPoint]=buffer_dis(Awake_logical,segmentLength,round(0.5*segmentLength),Electrical.points);
    [NREM_logseg,NREM_midPoint]  =buffer_dis(NREM_logical ,segmentLength,round(0.5*segmentLength),Electrical.points);
    [REM_logseg,REM_midPoint]    =buffer_dis(REM_logical  ,segmentLength,round(0.5*segmentLength),Electrical.points);
    logseg={Awake_logseg,NREM_logseg,REM_logseg};
    for ii=1:logseg_nState
        D6_Name{ii+(i-1)*logseg_nState}=E_Name;
        D6_Genotype{ii+(i-1)*logseg_nState}=E_Genotype;
        D6_Day{ii+(i-1)*logseg_nState}=E_Day;
        D6_State(ii+(i-1)*logseg_nState)=logseg_states(ii);
        
        logseg_nSegment=size(logseg{ii},2);
        D6_Data{ii+(i-1)*logseg_nState}=zeros(4,logseg_nSegment);
        D6_Data{ii+(i-1)*logseg_nState}(1,:)=mean(logseg{ii},1);
        for iii=1:logseg_nSegment
            logseg_points=logical2points(Awake_logseg(:,iii));
            logseg_points_lengths=(logseg_points(:,2)-logseg_points(:,1))/Electrical.fs;
            D6_Data{ii+(i-1)*logseg_nState}(2,iii)=mean(logseg_points_lengths);
            D6_Data{ii+(i-1)*logseg_nState}(3,iii)=std(logseg_points_lengths);
            D6_Data{ii+(i-1)*logseg_nState}(4,iii)=size(logseg_points,1);
        end
    end
    
    % DataTable2 with PSDs and length distributions
    temp_points={Electrical_points.Awake_Light;Electrical_points.NREM_Light;Electrical_points.REM_Light;...
                 Electrical_points.Awake_Dark;Electrical_points.NREM_Dark;Electrical_points.REM_Dark;};
    for ii=1:nState
        temp_points{ii}=temp_points{ii}(temp_points{ii}(:,2)-temp_points{ii}(:,1)+1>=FFT.nWindowLength,:);
    end
    temp_lengthdis={Awake_Light;NREM_Light;REM_Light;...
                    Awake_Dark;NREM_Dark;REM_Dark};
                
    % preallocation
%     PSD_cell=cell(nState*Settings.nChannel,1);
%     Name_cell=cell(nState*Settings.nChannel,1);
%     Genotype_cell=cell(nState*Settings.nChannel,1);
%     Day_cell=cell(nState*Settings.nChannel,1);
%     State_cell=cell(nState*Settings.nChannel,1);
%     Channel_cell=cell(nState*Settings.nChannel,1);
    clearvars PSD_cell Name_cell Genotype_cell Day_cell State_cell Channel_cell
    
%     Name_cell2=cell(nState,1);
%     Genotype_cell2=cell(nState,1);
%     Day_cell2=cell(nState,1);
%     State_cell2=cell(nState,1);
%     LengthDistribution_cell2=cell(nState,1);
    clearvars Name_cell2 Genotype_cell2 Day_cell2 State_cell2 LengthDistribution_cell2
                
    count1=0;
    count2=0;
    for ii=1:nState
        for iii=1:Settings.nChannel
            
            % for DataTable2
            if ~isempty(temp_points{ii})
            count1=count1+1;
            PSD_matrix=spectrogram_dis(Electrical.CH1234(:,iii),temp_points{ii},...
            FFT.nWindowLength,FFT.nOverlap,FFT.n,Electrical.fs);
            Name_cell{count1,1}=E_Name;
            Genotype_cell{count1,1}=E_Genotype;
            Day_cell{count1,1}=E_Day;
            State_cell{count1,1}=States{ii};
            Channel_cell{count1,1}=Settings.Channels{iii};
            PSD_cell{count1,1}=mean(PSD_matrix,2);
            end
            
        end
        
        % for DataTable3
        if ~isempty(temp_lengthdis{ii})
        count2=count2+1;
        Name_cell2{count2,1}=E_Name;
        Genotype_cell2{count2,1}=E_Genotype;
        Day_cell2{count2,1}=E_Day;
        State_cell2{count2,1}=States{ii};
        LengthDistribution_cell2{count2,1}=temp_lengthdis{ii};
        end
        
    end
    Name_Percentages_cell{i,1}=E_Name;
    Genotype_Percentages_cell{i,1}=E_Genotype;
    Day_Percentages_cell{i,1}=E_Day;
    Light_Percentages_cell{i,1}=Light_Percentages;
    Dark_Percentages_cell{i,1}=Dark_Percentages;
    
    DataTable2_cell{i,1}=table(Name_cell,Genotype_cell,Day_cell,State_cell,Channel_cell,PSD_cell,...
                    'VariableNames',{'Name','Genotype','Day','State','Channel','PSD'});
    DataTable3_cell{i,1}=table(Name_cell2,Genotype_cell2,Day_cell2,State_cell2,LengthDistribution_cell2,...
                    'VariableNames',{'Name','Genotype','Day','State','LengthDistribution'});
        
% States
States={'Light';'Dark'};
nState=length(States);
    
    % preallocation
%     PSD_cell=cell(nState*Settings.nChannel,1);
%     Name_cell=cell(nState*Settings.nChannel,1);
%     Genotype_cell=cell(nState*Settings.nChannel,1);
%     Day_cell=cell(nState*Settings.nChannel,1);
%     State_cell=cell(nState*Settings.nChannel,1);
%     Channel_cell=cell(nState*Settings.nChannel,1);
    clearvars PSD_cell Name_cell Genotype_cell Day_cell State_cell Channel_cell
    
    temp_points={points_Light;points_Dark};
    for ii=1:nState
        temp_points{ii}=temp_points{ii}(temp_points{ii}(:,2)-temp_points{ii}(:,1)+1>=FFT.nWindowLength,:);
    end
    
    count1=0;
    for ii=1:nState
        for iii=1:Settings.nChannel
            
            % for DataTable4
            if ~isempty(temp_points{ii})
            count1=count1+1;
            PSD_matrix=spectrogram_dis(Electrical.CH1234(:,iii),temp_points{ii},...
            FFT.nWindowLength,FFT.nOverlap,FFT.n,Electrical.fs);
            Name_cell{count1,1}=E_Name;
            Genotype_cell{count1,1}=E_Genotype;
            Day_cell{count1,1}=E_Day;
            State_cell{count1,1}=States{ii};
            Channel_cell{count1,1}=Settings.Channels{iii};
            PSD_cell{count1,1}=mean(PSD_matrix,2);
            end
            
        end
    end
    DataTable4_cell{i,1}=table(Name_cell,Genotype_cell,Day_cell,State_cell,Channel_cell,PSD_cell,...
                    'VariableNames',{'Name','Genotype','Day','State','Channel','PSD'});
                  
% States
States={'Awake';'NREM';'REM'};
nState=length(States);
    
    % preallocation
%     PSD_cell=cell(nState*Settings.nChannel,1);
%     Name_cell=cell(nState*Settings.nChannel,1);
%     Genotype_cell=cell(nState*Settings.nChannel,1);
%     Day_cell=cell(nState*Settings.nChannel,1);
%     State_cell=cell(nState*Settings.nChannel,1);
%     Channel_cell=cell(nState*Settings.nChannel,1);
    clearvars PSD_cell Name_cell Genotype_cell Day_cell State_cell Channel_cell
    
    % Create points
    points_Awake=logical2points((points2logical(Electrical.points_AA,Electrical.n)|points2logical(Electrical.points_QA,Electrical.n)));
    points_NREM=logical2points(points2logical(Electrical.points_NREM,Electrical.n));
    points_REM=logical2points(points2logical(Electrical.points_REM,Electrical.n));
    
    temp_points={points_Awake;points_NREM;points_REM};
    for ii=1:nState
        temp_points{ii}=temp_points{ii}(temp_points{ii}(:,2)-temp_points{ii}(:,1)+1>=FFT.nWindowLength,:);
    end
            
    count1=0;
    for ii=1:nState
        for iii=1:Settings.nChannel
            
            % for DataTable5
            if ~isempty(temp_points{ii})
            count1=count1+1;
            PSD_matrix=spectrogram_dis(Electrical.CH1234(:,iii),temp_points{ii},...
            FFT.nWindowLength,FFT.nOverlap,FFT.n,Electrical.fs);
            Name_cell{count1,1}=E_Name;
            Genotype_cell{count1,1}=E_Genotype;
            Day_cell{count1,1}=E_Day;
            State_cell{count1,1}=States{ii};
            Channel_cell{count1,1}=Settings.Channels{iii};
            PSD_cell{count1,1}=mean(PSD_matrix,2);
            end
            
        end
    end
    DataTable5_cell{i,1}=table(Name_cell,Genotype_cell,Day_cell,State_cell,Channel_cell,PSD_cell,...
                    'VariableNames',{'Name','Genotype','Day','State','Channel','PSD'});
    
end
close(bar)

DataTable_Percentages=table(Name_Percentages_cell,Genotype_Percentages_cell,Day_Percentages_cell,Light_Percentages_cell,Dark_Percentages_cell,...
                      'VariableNames',{'Name','Genotype','Day','LightPercentages','DarkPercentages'});
DataTable2=cat(1,DataTable2_cell{:});
DataTable3=cat(1,DataTable3_cell{:});
DataTable4=cat(1,DataTable4_cell{:});
DataTable5=cat(1,DataTable5_cell{:});
DataTable6=table(D6_Name,D6_Genotype,D6_Day,D6_State,D6_Data,...
           'VariableNames',{'Name','Genotype','Day','State','Data'});
save('DataTable','DataTable2','DataTable3','DataTable_Percentages','DataTable4','DataTable5','DataTable6')