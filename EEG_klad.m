% 1) Neurologger2 convert
% 2) MATLAB convert
% 3) Select
% 4) Analysis
% 5) Add starting time
% 6) Plot and statistics

%% number of genotypes in table
load('DataTable')
[~,b,~]=unique(DataTable.Name);
[a,~,c]=unique(DataTable(b,:).Genotype);
c=histcounts(c);
table(a,c')
%%
N=1000;
threshold=round(linspace(1,2500,N));
numberofrejectedseconds=zeros(N,4);
for i=1:N
numberofrejectedseconds(i,:)=sum(abs(Electrical.CH1234)>threshold(i));
end
numberofrejectedseconds=numberofrejectedseconds/Electrical.fs;
figure
plot(threshold,numberofrejectedseconds)

%%% OLD AVERAGE OVER AREAS
% Average over areas
count=0;
temptable.Area=strings(length(temptable.Channel),1);
[uniq.Channel,~,uniq.Channel_index]=unique(temptable.Channel);
nChannel=length(uniq.Channel);
for i=1:nChannel
    temp.string=strsplit(uniq.Channel{i},'-');
    temptable.Area(uniq.Channel_index==i)=temp.string(2);
end

uniq.Area=unique(temptable.Area);
nArea=length(uniq.Area);
for i=1:nArea
logical.area=strcmp(uniq.Area(i),temptable.Area);
uniq.State=unique(temptable.State(logical.area));
nState=length(uniq.State);
for ii=1:nState
logical.area_state=logical.area&strcmp(uniq.State(ii),temptable.State);
uniq.Name=unique(temptable.Name(logical.area_state));
nName=length(uniq.Name);
for iiii=1:nName
logical.area_state_name=logical.area_state&strcmp(uniq.Name(iiii),temptable.Name);

if sum(logical.area_state_name)
    
    % temporary 'table'
count=count+1;
bandtable.Area(count,1)=uniq.Area(i);
bandtable.State(count,1)=uniq.State(ii);
bandtable.Genotype(count,1)=uniq.Genotype(iii);
for bbb=1:Settings.nBand
bandtable.PSD{count,1}(bbb,:)=bandpower2(temp.power_arr,fff,Settings.BandRanges(bbb,:));
end
    
    
    count=count+1;
    temptable2.Area(count,1)=temptable.Area(find(logical.area_state_name,1,'first'));
    temptable2.State(count,1)=temptable.State(find(logical.area_state_name,1,'first'));
    temptable2.Genotype(count,1)=temptable.Genotype(find(logical.area_state_name,1,'first'));
    temptable2.PSD{count,1}=mean(cat(2,temptable.PSD{logical.area_state_name}),2);
end

end
end
end       
%% plot mean over right and left with meaned SEM (abs and norm)
% PSD/SEM_cell and area/State/temptable2.Genotype
uniq.State_full=unique(temptable2.State);
nState_full=length(uniq.State_full);
uniq.Genotype_full=unique(temptable2.Genotype);
Genotype_full=length(uniq.Genotype_full);

count=0;

uniq.Area=unique(temptable2.Area);
Area=length(uniq.Area);
for i=1:Area
    
    % 1 figure per area, save handles and legend names for each subplot in
    % each figure
    figure('Name',[uniq.Area{i},' Mean over Genotype'])
    legend_cell=cell(nState_full,2);
    handle_cell=cell(nState_full,2);
    
    logical.area=strcmp(temptable2.Area,uniq.Area(i));
    uniq.State=unique(temptable2.State(logical.area));
    nState=length(uniq.State);
    for ii=1:nState
        logical.area_state=logical.area&strcmp(temptable2.State,uniq.State(ii));
        uniq.Genotype=unique(temptable2.Genotype(logical.area_state));
        nGenotype=length(uniq.Genotype);
        for iii=1:nGenotype
            logical.area_state_genotype=logical.area_state&strcmp(temptable2.Genotype,uniq.Genotype(iii));
            
% start of content            
% state index
state_index=find(strcmp(uniq.State_full,uniq.State(ii)));
% genotype index
genotype_index=find(strcmp(uniq.Genotype_full,uniq.Genotype(iii)));

% absolute
temp.power_arr=cat(2,temptable2.PSD{logical.area_state_genotype});
[temp.power,temp.SEM]=meanSEM(temp.power_arr);
    
% absolute plot
subplot(nState_full,2,(state_index-1)*2+1)
h=plot(fff,temp.power,...
    'Color',c(genotype_index,:),'LineStyle',ls{genotype_index},'LineWidth',lw); hold on
patch([fff,...
       fff(end:-1:1),...
       fff(1)],...
      [temp.power-temp.SEM;...
       temp.power(end:-1:1)+temp.SEM(end:-1:1);...
       temp.power(1)],...
       c(genotype_index,:),...
       'FaceAlpha',SEM_transparency,...
       'EdgeAlpha',SEM_transparency)
title([uniq.Area{i},' ',uniq.State_full{state_index},' -absolute'])
xlim([Settings.HP_ele,maxFreq])
legend_cell{ii,1}=[legend_cell{ii,1},temptable2.Genotype(find(logical.area_state_genotype,1,'first'))];
handle_cell{ii,1}=[handle_cell{ii,1},h];
xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    
% normalization
temp.power_arr=(temp.power_arr./bandpower2(temp.power_arr,fff,[Settings.HP_ele,normFreq]))*100;
[temp.power_norm,temp.SEM_norm]=meanSEM(temp.power_arr);

% normalized plot
subplot(nState_full,2,state_index*2)
h=plot(fff,temp.power_norm,...
    'Color',c(genotype_index,:),'LineStyle',ls{genotype_index},'LineWidth',lw); hold on
patch([fff,...
       fff(end:-1:1),...
       fff(1)],...
      [temp.power_norm-temp.SEM_norm;...
       temp.power_norm(end:-1:1)+temp.SEM_norm(end:-1:1);...
       temp.power_norm(1)],...
       c(genotype_index,:),...
       'FaceAlpha',SEM_transparency,...
       'EdgeAlpha',SEM_transparency)
title([uniq.Area{i},' ',uniq.State_full{state_index},' -normalized'])
xlim([Settings.HP_ele,maxFreq])
legend_cell{ii,2}=[legend_cell{ii,2},temptable2.Genotype(find(logical.area_state_genotype,1,'first'))];
handle_cell{ii,2}=[handle_cell{ii,2},h];
xlabel('Frequency (Hz)'), ylabel('PSD (%)')
    
% temporary 'table'
count=count+1;
bandtable.Area(count,1)=uniq.Area(i);
bandtable.State(count,1)=uniq.State(ii);
bandtable.Genotype(count,1)=uniq.Genotype(iii);
for bbb=1:Settings.nBand
bandtable.PSD{count,1}(bbb,:)=bandpower2(temp.power_arr,fff,Settings.BandRanges(bbb,:));
end

        end
    end
    
    for l=1:nState_full
    for ll=1:2
        subplot(nState_full,2,ll+(l-1)*2)
        legend(handle_cell{l,ll},legend_cell{l,ll})
    end
    end
    
end


% %% Light/Dark assume to be from 7 AM to 7 PM, compute percentages and state length distributions
% % 'Table' strings
% E_Name=cell(nE,1);
% E_Genotype=cell(nE,1);
% E_Day=cell(nE,1);
% % Distributions
% Awake_Light=cell(nE,1);
% Sleep_Light=cell(nE,1);
% NREM_Light=cell(nE,1);
% REM_Light=cell(nE,1);
% Awake_Dark=cell(nE,1);
% Sleep_Dark=cell(nE,1);
% NREM_Dark=cell(nE,1);
% REM_Dark=cell(nE,1);
% % Percentages
% Light_Percentages=zeros(nE,4);
% Dark_Percentages=zeros(nE,4);
% 
% for i=1:nE
%     
%     % Load and create Light and Dark logicals
%     load(E_FileName{i},'header','Electrical')
%     starttimesplit=strsplit(header.starttime,'.');
%     StartTime=str2double(starttimesplit{1})*60*60+str2double(starttimesplit{2})*60+str2double(starttimesplit{3});
%     Electrical.t=mod(Electrical.t+StartTime,24*60*60);
%     Light_logical=Electrical.t>LightDark(1)*60*60&Electrical.t<LightDark(2)*60*60;
%     Dark_logical=~Light_logical;
%     
%     % Temptable strings
%     Esplit=strsplit(E_FileName{i},'_');
%     E_Name{i}=Esplit{1};
%     E_Genotype{i}=Esplit{2};
%     E_Day{i}=Esplit{3};
% 
%     % Distributions
%     % Light
%     Awake_Light_points=logical2points((points2logical(Electrical.points_AA,Electrical.n)|points2logical(Electrical.points_QA,Electrical.n))&Light_logical);
%     Awake_Light{i}=(Awake_Light_points(:,2)-Awake_Light_points(:,1)+1)/Electrical.fs;
%     Sleep_Light_points=logical2points((points2logical(Electrical.points_NREM,Electrical.n)|points2logical(Electrical.points_REM,Electrical.n))&Light_logical);
%     Sleep_Light{i}=(Sleep_Light_points(:,2)-Sleep_Light_points(:,1)+1)/Electrical.fs;
%     NREM_Light_points=logical2points(points2logical(Electrical.points_NREM,Electrical.n)&Light_logical);
%     NREM_Light{i}=(NREM_Light_points(:,2)-NREM_Light_points(:,1)+1)/Electrical.fs;
%     REM_Light_points=logical2points(points2logical(Electrical.points_REM,Electrical.n)&Light_logical);
%     REM_Light{i}=(REM_Light_points(:,2)-REM_Light_points(:,1)+1)/Electrical.fs;
%     % Dark
%     Awake_Dark_points=logical2points((points2logical(Electrical.points_AA,Electrical.n)|points2logical(Electrical.points_QA,Electrical.n))&Dark_logical);
%     Awake_Dark{i}=(Awake_Dark_points(:,2)-Awake_Dark_points(:,1)+1)/Electrical.fs;
%     Sleep_Dark_points=logical2points((points2logical(Electrical.points_NREM,Electrical.n)|points2logical(Electrical.points_REM,Electrical.n))&Dark_logical);
%     Sleep_Dark{i}=(Sleep_Dark_points(:,2)-Sleep_Dark_points(:,1)+1)/Electrical.fs;
%     NREM_Dark_points=logical2points(points2logical(Electrical.points_NREM,Electrical.n)&Dark_logical);
%     NREM_Dark{i}=(NREM_Dark_points(:,2)-NREM_Dark_points(:,1)+1)/Electrical.fs;
%     REM_Dark_points=logical2points(points2logical(Electrical.points_REM,Electrical.n)&Dark_logical);
%     REM_Dark{i}=(REM_Dark_points(:,2)-REM_Dark_points(:,1)+1)/Electrical.fs;
%     
%     % Percentages
%     Light_Percentages(i,:)=[sum(Awake_Light{i}),sum(Sleep_Light{i}),sum(NREM_Light{i}),sum(REM_Light{i})]/(sum(Light_logical)/Electrical.fs);
%     Dark_Percentages(i,:)=[sum(Awake_Dark{i}),sum(Sleep_Dark{i}),sum(NREM_Dark{i}),sum(REM_Dark{i})]/(sum(Dark_logical)/Electrical.fs);
%     
% end
% 
% %% Plot
% FIG.REM_Light=figure('Name','REM Light'); hold on
% FIG.REM_Dark=figure('Name','REM dark'); hold on
% 
% uniq.Name=unique(E_Name);
% nName=length(uniq.Name);
% REM_Light_Mean=zeros(nName,1);
% REM_Dark_Mean=zeros(nName,1);
% Genotype_Mean=cell(nName,1);
% for i=1:nName
%     logical.Name=strcmp(uniq.Name{i},E_Name);
%     set(0,'CurrentFigure',FIG.REM_Light)
%     temp=cat(1,REM_Light{logical.Name});
%     REM_Light_Mean(i)=mean(temp);
%     histogram(temp,'Normalization','pdf','BinWidth',5e3)
%     set(0,'CurrentFigure',FIG.REM_Dark)
%     temp=cat(1,REM_Dark{logical.Name});
%     REM_Dark_Mean(i)=mean(temp);
%     histogram(cat(1,REM_Dark{logical.Name}),'Normalization','pdf','BinWidth',5e3)
%     Genotype_Mean{i}=E_Genotype{find(logical.Name,1,'first')};
% end
% set(0,'CurrentFigure',FIG.REM_Light)
% legend(uniq.Name)
% set(0,'CurrentFigure',FIG.REM_Dark)
% legend(uniq.Name)
% 
% %% Plot Means
% uniq.Genotype=unique(Genotype_Mean);
% nGenotype=length(uniq.Genotype);
% FIG.REM_Light_Mean=figure('Name','REM Light Mean'); hold on
% FIG.REM_Dark_Mean=figure('Name','REM Dark Mean'); hold on
% for i=1:nGenotype
%     logical.Genotype=strcmp(uniq.Genotype{i},Genotype_Mean);
%     set(0,'CurrentFigure',FIG.REM_Light_Mean)
%     histogram(REM_Light_Mean(logical.Genotype))
%     set(0,'CurrentFigure',FIG.REM_Dark_Mean)
%     histogram(REM_Dark_Mean(logical.Genotype))
% end
% set(0,'CurrentFigure',FIG.REM_Light_Mean)
% legend(uniq.Genotype)
% set(0,'CurrentFigure',FIG.REM_Dark_Mean)
% legend(uniq.Genotype)
% 
% %%
% % for i=1:nGenotype
% % ttest2('vartype','unequal')


%% Check high CC
CC=zeros(6,nE);
for i=1:nE
    load(E(i).name,'Epoch')
    CC(:,i)=[mean(Epoch.CC{1,2}(:,end));...
             mean(Epoch.CC{1,3}(:,end));...
             mean(Epoch.CC{2,3}(:,end));...
             mean(Epoch.CC{1,4}(:,end));...
             mean(Epoch.CC{2,4}(:,end));...
             mean(Epoch.CC{3,4}(:,end))];
    if any(CC(:,i)>.99)
        load(E(i).name,'Electrical','FFT','Settings','Mouse')
        
        
temp_points={Electrical.points_AA,Electrical.points_NREM,Electrical.points_REM,Electrical.points_QA};
for iiii=1:Settings.nState
temp_points{iiii}=temp_points{iiii}(temp_points{iiii}(:,2)-temp_points{iiii}(:,1)+1>=FFT.nWindowLength,:);
end
PSD=cell(Settings.nChannel,Settings.nState);
for ii=1:Settings.nChannel
    figure
    for iiii=1:Settings.nState
        if ~isempty(temp_points{iiii})
        [PSD_matrix,t,f]=spectrogram_dis(Electrical.CH1234(:,ii),temp_points{iiii},...
         FFT.nWindowLength,FFT.nOverlap,FFT.n,Electrical.fs);
     subplot(4,1,iiii)
     imagesc(t,f,PSD_matrix)
     PSD{ii,iiii}=mean(PSD_matrix,2);
     title([Mouse,' ',Settings.Channels{ii},' ',Settings.State{iiii}],'Interpreter','none')
     xlabel('Time (s)','Interpreter','latex'),ylabel('Frequency (Hz)','Interpreter','latex')
     set(gca,'YDir','normal')
     set(gca,'CLim',[0,mean(max(PSD_matrix))])
     ylim([0,20])
        end
    end
end

for ii=1:Settings.nChannel
count=0;
    figure
    hold on
    grid on
%     plot_h=gobjects(Settings.nChannel,1);
    for iiii=1:Settings.nState
        if ~isempty(temp_points{iiii})
            count=count+1;
    plot_h(count)=plot(f,PSD{ii,iiii},'Color',Settings.c(iiii,:),'LineWidth',2);
    legend_cell{count}=Settings.State{iiii};
        end
    end
    title([Mouse,' ',Settings.Channels{ii}],'Interpreter','none')
    xlabel('Frequency (Hz)','Interpreter','latex'),ylabel('Power ($\mu V^2$)','Interpreter','latex')
    xlim([0,50])
    legend(plot_h,legend_cell)
end
        
        
        
    end
end



%% OLD EEG plot
% Plot mean PSD for each mouse and genotype
clearvars, close all
%% Import
EF = uigetdir('','Select Experiment Folder (EF)');
cd(EF)
E=dir('*_Analyzed.mat');
nE=length(E);
DataTable_cell1=cell(nE,1);
DataTable_cell2=cell(nE,1);
DataTable_cell3=cell(nE,1);
E_FileName=cell(nE,1);
for i=1:nE
    E_FileName{i}=E(i).name(1:end-4);
    load(E_FileName{i},'DataTable','DataTable2','DataTable3')
    nM=size(DataTable,1);
    checkifempty=true(nM,1);
    for ii=1:nM
        checkifempty(ii)=~isempty(DataTable.Power{ii});
    end
    
    DataTable_cell1{i}=DataTable(checkifempty,:);
    DataTable_cell2{i}=DataTable2;
    DataTable_cell3{i}=DataTable3;
end
DataTable1=cat(1,DataTable_cell1{:});
DataTable2=cat(1,DataTable_cell2{:});
DataTable3=cat(1,DataTable_cell3{:});

load('Settings')

% Leave out high cross correlation between channels


%% plot mean over right and left with meaned SEM (abs and norm)
% Power/SEM_cell and area/State/temptable2.Genotype
uniq.State_full=unique(temptable2.State);
nState_full=length(uniq.State_full);
uniq.Genotype_full=unique(temptable2.Genotype);
Genotype_full=length(uniq.Genotype_full);

count=0;

uniq.Area=unique(temptable2.Area);
Area=length(uniq.Area);
for i=1:Area
    
    % 1 figure per area, save handles and legend names for each subplot in
    % each figure
    figure('Name',[uniq.Area{i},' Mean over Genotype'])
    legend_cell=cell(nState_full,2);
    handle_cell=cell(nState_full,2);
    
    logical.area=strcmp(temptable2.Area,uniq.Area(i));
    uniq.State=unique(temptable2.State(logical.area));
    nState=length(uniq.State);
    for ii=1:nState
        logical.area_state=logical.area&strcmp(temptable2.State,uniq.State(ii));
        uniq.Genotype=unique(temptable2.Genotype(logical.area_state));
        nGenotype=length(uniq.Genotype);
        for iii=1:nGenotype
            logical.area_state_genotype=logical.area_state&strcmp(temptable2.Genotype,uniq.Genotype(iii));
            
% start of content            
% state index
state_index=find(strcmp(uniq.State_full,uniq.State(ii)));
% genotype index
genotype_index=find(strcmp(uniq.Genotype_full,uniq.Genotype(iii)));

% absolute
temp.power_arr=cat(2,temptable2.Power{logical.area_state_genotype});
[temp.power,temp.SEM]=meanSEM(temp.power_arr);
    
% absolute plot
subplot(nState_full,2,(state_index-1)*2+1)
h=plot(fff,temp.power,...
    'Color',c(genotype_index,:),'LineStyle',ls{genotype_index},'LineWidth',lw); hold on
patch([fff,...
       fff(end:-1:1),...
       fff(1)],...
      [temp.power-temp.SEM;...
       temp.power(end:-1:1)+temp.SEM(end:-1:1);...
       temp.power(1)],...
       c(genotype_index,:),...
       'FaceAlpha',SEM_transparency,...
       'EdgeAlpha',SEM_transparency)
title([uniq.Area{i},' ',uniq.State_full{state_index},' -absolute'])
xlim([Settings.HP_ele,maxFreq])
legend_cell{ii,1}=[legend_cell{ii,1},temptable2.Genotype(find(logical.area_state_genotype,1,'first'))];
handle_cell{ii,1}=[handle_cell{ii,1},h];
xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    
% normalization
temp.power_arr=(temp.power_arr./bandpower2(temp.power_arr,fff,[Settings.HP_ele,normFreq]))*100;
[temp.power_norm,temp.SEM_norm]=meanSEM(temp.power_arr);

% normalized plot
subplot(nState_full,2,state_index*2)
h=plot(fff,temp.power_norm,...
    'Color',c(genotype_index,:),'LineStyle',ls{genotype_index},'LineWidth',lw); hold on
patch([fff,...
       fff(end:-1:1),...
       fff(1)],...
      [temp.power_norm-temp.SEM_norm;...
       temp.power_norm(end:-1:1)+temp.SEM_norm(end:-1:1);...
       temp.power_norm(1)],...
       c(genotype_index,:),...
       'FaceAlpha',SEM_transparency,...
       'EdgeAlpha',SEM_transparency)
title([uniq.Area{i},' ',uniq.State_full{state_index},' -normalized'])
xlim([Settings.HP_ele,maxFreq])
legend_cell{ii,2}=[legend_cell{ii,2},temptable2.Genotype(find(logical.area_state_genotype,1,'first'))];
handle_cell{ii,2}=[handle_cell{ii,2},h];
xlabel('Frequency (Hz)'), ylabel('PSD (%)')
    
% temporary 'table'
count=count+1;
bandtable.Area(count,1)=uniq.Area(i);
bandtable.State(count,1)=uniq.State(ii);
bandtable.Genotype(count,1)=uniq.Genotype(iii);
for bbb=1:Settings.nBand
bandtable.Power{count,1}(bbb,:)=bandpower2(temp.power_arr,fff,Settings.BandRanges(bbb,:));
end

        end
    end
    
    for l=1:nState_full
    for ll=1:2
        subplot(nState_full,2,ll+(l-1)*2)
        legend(handle_cell{l,ll},legend_cell{l,ll})
    end
    end
    
end

%% Light/Dark assume to be from 7 AM to 7 PM, compute percentages and state length distributions
% 'Table' strings
E_Name=cell(nE,1);
E_Genotype=cell(nE,1);
E_Day=cell(nE,1);
% Distributions
Awake_Light=cell(nE,1);
Sleep_Light=cell(nE,1);
NREM_Light=cell(nE,1);
REM_Light=cell(nE,1);
Awake_Dark=cell(nE,1);
Sleep_Dark=cell(nE,1);
NREM_Dark=cell(nE,1);
REM_Dark=cell(nE,1);
% Percentages
Light_Percentages=zeros(nE,4);
Dark_Percentages=zeros(nE,4);

for i=1:nE
    
    % Load and create Light and Dark logicals
    load(E_FileName{i},'header','Electrical')
    starttimesplit=strsplit(header.starttime,'.');
    StartTime=str2double(starttimesplit{1})*60*60+str2double(starttimesplit{2})*60+str2double(starttimesplit{3});
    Electrical.t=mod(Electrical.t+StartTime,24*60*60);
    Light_logical=Electrical.t>LightDark(1)*60*60&Electrical.t<LightDark(2)*60*60;
    Dark_logical=~Light_logical;
    
    % Temptable strings
    Esplit=strsplit(E_FileName{i},'_');
    E_Name{i}=Esplit{1};
    E_Genotype{i}=Esplit{2};
    E_Day{i}=Esplit{3};

    % Distributions
    % Light
    Awake_Light_points=logical2points((points2logical(Electrical.points_AA,Electrical.n)|points2logical(Electrical.points_QA,Electrical.n))&Light_logical);
    Awake_Light{i}=(Awake_Light_points(:,2)-Awake_Light_points(:,1)+1)/Electrical.fs;
    Sleep_Light_points=logical2points((points2logical(Electrical.points_NREM,Electrical.n)|points2logical(Electrical.points_REM,Electrical.n))&Light_logical);
    Sleep_Light{i}=(Sleep_Light_points(:,2)-Sleep_Light_points(:,1)+1)/Electrical.fs;
    NREM_Light_points=logical2points(points2logical(Electrical.points_NREM,Electrical.n)&Light_logical);
    NREM_Light{i}=(NREM_Light_points(:,2)-NREM_Light_points(:,1)+1)/Electrical.fs;
    REM_Light_points=logical2points(points2logical(Electrical.points_REM,Electrical.n)&Light_logical);
    REM_Light{i}=(REM_Light_points(:,2)-REM_Light_points(:,1)+1)/Electrical.fs;
    % Dark
    Awake_Dark_points=logical2points((points2logical(Electrical.points_AA,Electrical.n)|points2logical(Electrical.points_QA,Electrical.n))&Dark_logical);
    Awake_Dark{i}=(Awake_Dark_points(:,2)-Awake_Dark_points(:,1)+1)/Electrical.fs;
    Sleep_Dark_points=logical2points((points2logical(Electrical.points_NREM,Electrical.n)|points2logical(Electrical.points_REM,Electrical.n))&Dark_logical);
    Sleep_Dark{i}=(Sleep_Dark_points(:,2)-Sleep_Dark_points(:,1)+1)/Electrical.fs;
    NREM_Dark_points=logical2points(points2logical(Electrical.points_NREM,Electrical.n)&Dark_logical);
    NREM_Dark{i}=(NREM_Dark_points(:,2)-NREM_Dark_points(:,1)+1)/Electrical.fs;
    REM_Dark_points=logical2points(points2logical(Electrical.points_REM,Electrical.n)&Dark_logical);
    REM_Dark{i}=(REM_Dark_points(:,2)-REM_Dark_points(:,1)+1)/Electrical.fs;
    
    % Percentages
    Light_Percentages(i,:)=[sum(Awake_Light{i}),sum(Sleep_Light{i}),sum(NREM_Light{i}),sum(REM_Light{i})]/(sum(Light_logical)/Electrical.fs);
    Dark_Percentages(i,:)=[sum(Awake_Dark{i}),sum(Sleep_Dark{i}),sum(NREM_Dark{i}),sum(REM_Dark{i})]/(sum(Dark_logical)/Electrical.fs);
    
end



%% 
E=dir('C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Analyzed\*_*_*_Analyzed.mat');
cd('C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Analyzed')
nE=length(E);
CC=zeros(6,nE);
for i=1:nE
    load(E(i).name,'Epoch')
    CC(:,i)=[mean(Epoch.CC{1,2}(:,end));...
             mean(Epoch.CC{1,3}(:,end));...
             mean(Epoch.CC{2,3}(:,end));...
             mean(Epoch.CC{1,4}(:,end));...
             mean(Epoch.CC{2,4}(:,end));...
             mean(Epoch.CC{3,4}(:,end))];
end
CC>.99
{E(any(CC>.99)).name}'

% 
% 
% 
% 
% 
% %% PRELIM fit 3Dn3
% 
% X=[Epoch.DynAcc,Epoch.CC{1,4}(:,2),Epoch.BandPowers{2}(:,6)];
% 
% GMM=fitgmdist(X,3);
% index=cluster(GMM,X);
% C=zeros(length(index),3);
% C(:,1)=index==1;
% C(:,2)=index==2;
% C(:,3)=index==3;
% 
% figure
% scatter3(X(:,1),X(:,2),X(:,3),5,C)
% xlabel('Length Acceleration (g)','Interpreter','latex')
% ylabel('$\theta$ L/R S CC','Interpreter','latex')
% zlabel('$\gamma$ power ($\mu V^2$)','Interpreter','latex')
% [sum(index==1),sum(index==2),sum(index==3)]/sum([sum(index==1),sum(index==2),sum(index==3)])
% 
% 
% 
% 
% 
% %% Theta cross 
% 
% GMMthetacross=fitgmdist(Epoch.CC{1,4}(:,2),2,'Options',Settings.Options,'Replicates',50);
% thetacross=linspace(0,1,1e5)';
% pdfthetacross=pdf(GMMthetacross,thetacross);
% [~,~,pthetacross_nonindexed]=cluster(GMMthetacross,Epoch.CC{1,4}(:,2));
% [~,INDEX]=sort(GMMthetacross.mu,'descend');
% INDEX_high=INDEX(1);
% INDEX_low=INDEX(2);
% pthetacross=pthetacross_nonindexed(:,[INDEX_low,INDEX_high]);
% clearvars pthetacross_nonindexed
%        
% COLOR  =   [0         0.4470    0.7410;
%             0.8500    0.3250    0.0980];
% MARKER =   {'o';'x'};
% 
% figure,histogram(Epoch.CC{1,4}(:,2),'Normalization','pdf')
% for i=1:Epoch.n
%     C(i,:)=pthetacross(i,:)*COLOR;
% end
% hold on
% grid on
% plot(thetacross,pdfthetacross,'LineWidth',2)
% title('\theta cross-correlation between left and right somatosensory cortex has mixed Gaussian distribution','Interpreter','tex')
% xlabel('\theta cross-correlation')
% ylabel('Epoch count (%)')
% xlim([.8,1])
% 
% for i=1:Settings.nChannel
% figure
% hold on
% grid on
% [~,index]=max(pthetacross,[],2);
% thetalow=index==1;
% thetahigh=index==2;
% State={thetalow,thetahigh};
% nState=length(State);
% Dec3=gobjects(nState,1);
% for iiii=1:nState
% Dec3(iiii)=scatter3(Epoch.DynAcc(State{iiii}),Epoch.DELTA{i}(State{iiii}),Epoch.THETA{i}(State{iiii}),5,C(State{iiii},:),...
% 'Marker',MARKER{iiii},'LineWidth',1);
% end
% title(Settings.Channels{i},'Interpreter','latex')
% xlabel('Length Acceleration (g)','Interpreter','latex')
% ylabel('NREM ratio ($\Delta=\frac{\delta*\alpha}{\beta*\gamma}$)','Interpreter','latex')
% zlabel('REM ratio ($\Theta=\frac{\theta^2}{\delta*\alpha}$)','Interpreter','latex')
% legend(Dec3,{['Low \Theta cross-correlation, mean=',num2str(GMMthetacross.mu(1))];['High \Theta cross-correlation, mean=',num2str(GMMthetacross.mu(2))]})
% set(gca,'CameraPosition',[-3.0596  -66.8125    8.2449])
% end
% 
% 
% 
% 
% 
% 
% 
% %% Combine overlapping probabilities in probability state x time matrix for left-M and right-M
% P=cell(Settings.nChannel,1);
% P3D=cell(Settings.nChannel,1);
% for i=1:Settings.nChannel
%     P{i}=zeros(Electrical.n,Settings.nState);
%     P3D{i}=zeros(Electrical.n,Settings.nState);
%     for iiii=1:Settings.nState
%         P{i}(:,iiii)=buffer_inv_dis(repmat(Epoch.p{i}(:,iiii)',[Electrical.nEpochLength,1]),...
%             Electrical.nEpochOverlap,'mean',Electrical.points,Electrical.n,NaN);
%         P3D{i}(:,iiii)=buffer_inv_dis(repmat(Epoch.p_3Dn4{i}(:,iiii)',[Electrical.nEpochLength,1]),...
%             Electrical.nEpochOverlap,'mean',Electrical.points,Electrical.n,NaN);
%     end
%     P{i}=P{i}./sum(P{i},2);
%     P3D{i}=P3D{i}./sum(P3D{i},2);
% end
% 
% %%% 3DModel end
% 
% % %% Fit 4 3D GM Models 
% % Epoch.p_3Dn4=cell(Settings.nChannel,1);
% % Epoch.ind_p_3Dn4=cell(Settings.nChannel,1);
% % GMModel_3Dn4=cell(Settings.nChannel,1);
% % Epoch.AA_3Dn4     = cell(Settings.nChannel,1);
% % Epoch.nAA_3Dn4    = zeros(Settings.nChannel,1);
% % Epoch.NREM_3Dn4   = cell(Settings.nChannel,1);
% % Epoch.nNREM_3Dn4  = zeros(Settings.nChannel,1);
% % Epoch.REM_3Dn4    = cell(Settings.nChannel,1);
% % Epoch.nREM_3Dn4   = zeros(Settings.nChannel,1);
% % Epoch.QA_3Dn4     = cell(Settings.nChannel,1);
% % Epoch.nQA_3Dn4    = zeros(Settings.nChannel,1);
% % for i=1:Settings.nChannel
% %     
% % [~,Suggestion]=max(Epoch.p{i},[],2);
% % X=[Epoch.DynAcc,Epoch.DELTA{i},Epoch.THETA{i}];
% % Xdup=cell(Settings.nChannel,1);
% % Suggestiondup=cell(Settings.nChannel,1);
% % for iiii=1:Settings.nState
% %     Xdup{iiii}=DuplicatePointsWithProbability(X,Epoch.p{i}(:,iiii),Settings.nMaxPoint);
% %     Suggestiondup{iiii}=DuplicatePointsWithProbability(Suggestion,Epoch.p{i}(:,iiii),Settings.nMaxPoint);
% % end
% % Xdup=cat(1,Xdup{:});
% % Suggestiondup=cat(1,Suggestiondup{:});
% % 
% % GMModel3Dn4_cell=cell(Settings.nChannel,1);
% % for iiii=1:Settings.nState
% % GMModel3Dn4_cell{iiii}=fitgmdist(Xdup(Suggestiondup==iiii,:),1,'Options',Settings.Options);
% % end
% % GMModel_3Dn4{i}=gmdistribution(cat(1,GMModel3Dn4_cell{1}.mu,GMModel3Dn4_cell{2}.mu,GMModel3Dn4_cell{3}.mu,GMModel3Dn4_cell{4}.mu),...
% %     cat(3,GMModel3Dn4_cell{1}.Sigma,GMModel3Dn4_cell{2}.Sigma,GMModel3Dn4_cell{3}.Sigma,GMModel3Dn4_cell{4}.Sigma),...
% %     cat(1,GMModel3Dn4_cell{1}.ComponentProportion,GMModel3Dn4_cell{2}.ComponentProportion,GMModel3Dn4_cell{3}.ComponentProportion,GMModel3Dn4_cell{4}.ComponentProportion));
% % [~,~,Epoch.p_3Dn4{i}]=cluster(GMModel_3Dn4{i},X);
% % [~,Epoch.ind_p_3Dn4{i}]=max(Epoch.p_3Dn4{i},[],2);
% % Epoch.AA_3Dn4{i}     = Epoch.ind_p_3Dn4{i}==1;
% % Epoch.nAA_3Dn4(i)    = sum(Epoch.AA_3Dn4{i});
% % Epoch.NREM_3Dn4{i}   = Epoch.ind_p_3Dn4{i}==2;
% % Epoch.nNREM_3Dn4(i)  = sum(Epoch.NREM_3Dn4{i});
% % Epoch.REM_3Dn4{i}    = Epoch.ind_p_3Dn4{i}==3;
% % Epoch.nREM_3Dn4(i)   = sum(Epoch.REM_3Dn4{i});
% % Epoch.QA_3Dn4{i}     = Epoch.ind_p_3Dn4{i}==4;
% % Epoch.nQA_3Dn4(i)    = sum(Epoch.QA_3Dn4{i});
% % 
% % for iii=1:nEpoch
% %     Epoch.C{i}(iii,:)=Epoch.p_3Dn4{i}(iii,:)*Settings.c(1:Settings.nState,:);
% % end
% % 
% % figure
% % hold on
% % grid on
% % Final=gobjects(Settings.nState,1);
% % State={Epoch.AA_3Dn4;Epoch.NREM_3Dn4;Epoch.REM_3Dn4;Epoch.QA_3Dn4};
% % for iiii=1:Settings.nState
% % Final(iiii)=scatter3(Epoch.DynAcc(State{iiii}{i}),Epoch.DELTA{i}(State{iiii}{i}),Epoch.THETA{i}(State{iiii}{i}),Epoch.S(State{iiii}{i}),Epoch.C{i}(State{iiii}{i},:),...
% % 'Marker',Settings.Marker{iiii},'LineWidth',1);
% % end
% % title(Settings.Channels{i})
% % xlabel('Length of Dynamic Acceleration vector (g)')
% % ylabel('NREM ratio ($\Delta=\frac{\delta*\alpha}{\beta*\gamma}$)')
% % zlabel('REM ratio ($\Theta=\frac{\theta^2}{\delta*\alpha}$)')
% % legend(Final,{'AA';'NREM';'REM';'QA'})
% % 
% % end
% 
% %%% 3DModel end
% 
% P3D_FA_ele=P3D{2}.*P3D{3};
% P3D_FA_ele=P3D_FA_ele./sum(P3D_FA_ele,2);
% 
% 
% Electrical.points3D_FA=logical2points(~isnan(P3D_FA_ele(:,1)));
% [~,P3D_FA_index]=max(P3D_FA_ele,[],2);
% Electrical.points3D_AA=logical2points(P3D_FA_index==1);
% Electrical.points3D_NREM=logical2points(P3D_FA_index==2);
% Electrical.points3D_REM=logical2points(P3D_FA_index==3);
% Electrical.points3D_QA=logical2points(P3D_FA_index==4);
% clearvars P3D_FA_index
% 
% 
% Electrical.points3D_FA(:,1)   =Electrical.points3D_FA(:,1)-1; Electrical.points3D_FA=round(Electrical.points3D_FA/EleAccMultiple)*EleAccMultiple; Electrical.points3D_FA(:,1)=Electrical.points3D_FA(:,1)+1;
% Electrical.points3D_AA(:,1)   =Electrical.points3D_AA(:,1)-1; Electrical.points3D_AA=round(Electrical.points3D_AA/EleAccMultiple)*EleAccMultiple; Electrical.points3D_AA(:,1)=Electrical.points3D_AA(:,1)+1;
% Electrical.points3D_NREM(:,1) =Electrical.points3D_NREM(:,1)-1; Electrical.points3D_NREM=round(Electrical.points3D_NREM/EleAccMultiple)*EleAccMultiple; Electrical.points3D_NREM(:,1)=Electrical.points3D_NREM(:,1)+1;
% Electrical.points3D_REM(:,1)  =Electrical.points3D_REM(:,1)-1; Electrical.points3D_REM=round(Electrical.points3D_REM/EleAccMultiple)*EleAccMultiple; Electrical.points3D_REM(:,1)=Electrical.points3D_REM(:,1)+1;
% Electrical.points3D_QA(:,1)   =Electrical.points3D_QA(:,1)-1; Electrical.points3D_QA=round(Electrical.points3D_QA/EleAccMultiple)*EleAccMultiple; Electrical.points3D_QA(:,1)=Electrical.points3D_QA(:,1)+1;
% 
% 
% Acceleration.points3D_FA   =[Electrical.points3D_FA(:,1)+EleAccMultiple-1,Electrical.points3D_FA(:,2)]/EleAccMultiple;
% Acceleration.points3D_AA   =[Electrical.points3D_AA(:,1)+EleAccMultiple-1,Electrical.points3D_AA(:,2)]/EleAccMultiple;
% Acceleration.points3D_NREM =[Electrical.points3D_NREM(:,1)+EleAccMultiple-1,Electrical.points3D_NREM(:,2)]/EleAccMultiple;
% Acceleration.points3D_REM  =[Electrical.points3D_REM(:,1)+EleAccMultiple-1,Electrical.points3D_REM(:,2)]/EleAccMultiple;
% Acceleration.points3D_QA   =[Electrical.points3D_QA(:,1)+EleAccMultiple-1,Electrical.points3D_QA(:,2)]/EleAccMultiple;
% 
% % %% Try to do 3D n=4 gaussian fit
% % X=[Epoch.DynAcc,Epoch.DELTA{3},Epoch.THETA{3}];
% % GMModel3Dn4=fitgmdist(X,4);
% % idx=cluster(GMModel3Dn4,X);
% % figure
% % scatter3(X(:,1),X(:,2),X(:,3),Epoch.S,Settings.c(idx,:))
% % 
% % %% Fit 4 final 3D GMModel with n=4 clusters by using previously calculated state as suggestions
% % Epoch.p_3Dn4=cell(nChannel,1);
% % Epoch.ind_p_3Dn4=cell(nChannel,1);
% % for i=1:nChannel
% %     
% % [~,Suggestion]=max(Epoch.p{i},[],2);
% % X=[Epoch.DynAcc,Epoch.DELTA{i},Epoch.THETA{i}];
% % Xdup=cell(nChannel,1);
% % Suggestiondup=cell(nChannel,1);
% % for iiii=1:nState_Full
% %     Xdup{iiii}=DuplicatePointsWithProbability(X,Epoch.p{i}(:,iiii),Settings.nMaxPoint);
% %     Suggestiondup{iiii}=DuplicatePointsWithProbability(Suggestion,Epoch.p{i}(:,iiii),Settings.nMaxPoint);
% % end
% % Xdup=cat(1,Xdup{:});
% % Suggestiondup=cat(1,Suggestiondup{:});
% % 
% % GMModel3Dn4=fitgmdist(Xdup,4,'Options',Settings.Options,'Start',Suggestiondup);
% % [~,~,Epoch.p_3Dn4{i}]=cluster(GMModel3Dn4,X);
% % [~,Epoch.ind_p_3Dn4{i}]=max(Epoch.p_3Dn4{i},[],2);
% % Epoch.AA{i}     = Epoch.ind_REM_QA{i}==1;
% % Epoch.nAA(i)    = sum(Epoch.AA{i});
% % Epoch.NREM{i}   = Epoch.ind_REM_QA{i}==2;
% % Epoch.nNREM(i)  = sum(Epoch.NREM{i});
% % Epoch.REM{i}    = Epoch.ind_REM_QA{i}==3;
% % Epoch.nREM(i)   = sum(Epoch.REM{i});
% % Epoch.QA{i}     = Epoch.ind_REM_QA{i}==4;
% % Epoch.nQA(i)    = sum(Epoch.QA{i});
% % 
% % for iii=1:nEpoch
% %     Epoch.C{i}(iii,:)=Epoch.p_3Dn4{i}(iii,:)*Settings.c(1:nState_Full,:);
% % end
% % 
% % figure
% % hold on
% % grid on
% % Final=gobjects(nState_Full,1);
% % State={Epoch.AA;Epoch.NREM;Epoch.REM;Epoch.QA};
% % for iiii=1:nState_Full
% % Final(iiii)=scatter3(Epoch.DynAcc(State{iiii}{i}),Epoch.DELTA{i}(State{iiii}{i}),Epoch.THETA{i}(State{iiii}{i}),Epoch.S(State{iiii}{i}),Epoch.C{i}(State{iiii}{i},:),...
% %     'Marker',Settings.Marker{iiii},'LineWidth',1);
% % end
% % title(Settings.Channels{i})
% % xlabel('Length of Dynamic Acceleration vector (g)')
% % ylabel('NREM ratio ($\Delta=\frac{\delta*\alpha}{\beta*\gamma}$)')
% % zlabel('REM ratio ($\Theta=\frac{\theta^2}{\delta*\alpha}$)')
% % legend(Final,{'AA';'NREM';'REM';'QA'})
% % 
% % end
% 
% % %% Altered final plot
% % % preallocation
% % Epoch.Marker=cell(nChannel,1);
% % for i=1:nChannel
% % Epoch.Marker{i}=cell(nEpoch,1);
% % end
% % 
% % for i=1:nChannel
% %     
% % figure
% % hold on
% % grid on
% % 
% % Epoch.Marker{i}(Epoch.AA{i})={'x'};
% % Epoch.Marker{i}(Epoch.NREM)={'o'};
% % Epoch.Marker{i}(Epoch.REM)={'*'};
% % Epoch.Marker{i}(Epoch.QA)={'^'};
% % 
% % for iii=1:nEpoch
% %     plot3(Epoch.DynAcc(iii),Epoch.DELTA{i}(iii),Epoch.THETA{i}(iii),'Color',Epoch.p_4{i}(iii,1:3),'Marker',Epoch.Marker{i}{iii})
% % end
% % 
% % h(1)=plot(NaN,'Color',[1,0,0]);
% % h(2)=plot(NaN,'Color',[0,1,0]);
% % h(3)=plot(NaN,'Color',[0,0,1]);
% % legend(h,{'AA','NREM','REM'})
% % 
% % title(Settings.Channels{i})
% % xlabel('Length of Dynamic Acceleration vector (g)')
% % ylabel('\Delta')
% % zlabel('\Theta')
% % 
% % end
% 
% 
% 
% % Epoch.logical_2=Epoch.p_2>1-Settings.p_value;
% 
% for iiii=1:CC.nWindow
%                 [pks,locs]=findpeaks(cc(:,iiii),'WidthReference','halfheight');
%                 if ~isempty(pks)
%                     [~,closestpeakto0_index]=min(abs(locs));
%                     CC_dis(iiii)  = pks(closestpeakto0_index);
%                     Lag_dis(iiii) = locs(closestpeakto0_index);
%                 else
%                     CC_dis(iiii)=cc(CC.nMaxLag+1,iiii);
%                     Lag_dis(iiii)=0;
%                 end
%             end
% 
% %% Compute overall RMS, CC and delay for each band and channel (3*nChannel*nBand dimensions) (for checking purposes)
% BandOverallRanges=[Settings.BandRanges(1,1),Settings.BandRanges(end,2);Settings.HP_ele,Settings.LP_ele];
% nBandOverall=size(BandOverallRanges,1);
% Epoch.RMS_Overall=cell(nChannel,1); % preallocation
% Epoch.CC_Overall=cell(nChannel,1); % preallocation
% Epoch.Delay_Overall=cell(nChannel,1); % preallocation
% for i=1:nChannel
%     Epoch.RMS_Overall{i}=zeros(nEpoch,nBandOverall)
%     for ii=1:nBandOverall
% Electrical.Epoch = HighLowPassfilter(Settings.order_ele, BandOverallRanges, Electrical.fs, Electrical.CH1234(:,i));
% Electrical.Epoch = buffer_dis(Electrical.Epoch,Electrical.nEpochLength,Electrical.nEpochOverlap,Electrical.points);
% Epoch.RMS{i}(:,ii) = rms(Electrical.Epoch);
%     end
% end
% 
% %% Band indices
% for i=1:nBand
%     eval([Settings.Bands{i},'=i;'])
% end
% 
% %% Plot
% figure
% X=[Bin.a,Bin.BandPowers{1}(:,1),Bin.BandPowers{1}(:,2)];
% scatter3(X(:,1),X(:,2),X(:,3))
% xlabel('a'),ylabel('\delta'),zlabel('\theta')
% 
% GMModel = fitgmdist(X(:,1),2);
% idx1=cluster(GMModel,X(:,1));
% % idx1=kmeans(X,2);
% 
% S=ones(nBin,1);
% C=Settings.c(idx1,:);
% figure
% scatter3(X(:,1),X(:,2),X(:,3),S,C)
% 
% [~,index]=sort(GMModel.mu,'descend');
% index_Awake=index(1);
% index_Sleep=index(2);
% logical_Sleep=idx1==index_Sleep;
% GMModel2 = fitgmdist(X(logical_Sleep),2);
% idx2=cluster(GMModel,X(logical_Sleep));
% 
% %%%%%%%%%%%%%%%%%%%%
% 
% Filtering of Accelerometer data and Truncation
% % create time vectors (to also truncate)
% Electrical.t=0:1/Electrical.fs:(size(Electrical.CH1234,1)-1)/Electrical.fs; % s
% Acceleration.t=0:1/Acceleration.fs:(size(Acceleration.XYZ,1)-1)/Acceleration.fs; % s
% 
% % truncate acceleration data with logical made by Select.m
% Acceleration.XYZ=Acceleration.XYZ(DataTable_Logical.Select_Acceleration{temp.select_index},:);
% Acceleration.t=Acceleration.t(DataTable_Logical.Select_Acceleration{temp.select_index});
% % truncate electrical data
% Electrical.CH1234=Electrical.CH1234(DataTable_Logical.Select_Electrical{temp.select_index},:);
% Electrical.t=Electrical.t(DataTable_Logical.Select_Electrical{temp.select_index});
% 
% % additional small truncation with modulus so that we have whole time bins
% % (e.g. 20-60 s)
% % truncate acceleration data
% Acceleration.nBin=round(Settings.timeBin*Acceleration.fs);
% temp.mod=mod(length(Acceleration.t),Acceleration.nBin);
% Acceleration.XYZ=Acceleration.XYZ(1:end-temp.mod,:);
% Acceleration.t=Acceleration.t(1:end-temp.mod);
% % truncate electrical data
% Electrical.nBin=round(Settings.timeBin*Electrical.fs);
% temp.mod=mod(length(Electrical.t),Electrical.nBin);
% Electrical.CH1234=Electrical.CH1234(1:end-temp.mod,:);
% Electrical.t=Electrical.t(1:end-temp.mod);
% 
% % some constants for acceleration data
% Acceleration.n=length(Acceleration.t);
% Acceleration.nCol=Acceleration.n/Acceleration.nBin;
% Acceleration.nDim=size(Acceleration.XYZ,2);
% % some constants for electrical data
% Electrical.n=length(Electrical.t);
% Electrical.nCol=Electrical.n/Electrical.nBin;
% Electrical.nChannel=size(Electrical.CH1234,2);
% 
% 
% %     % spec constants
% %     temp.fff=0:Acceleration.fs/temp.nfft:Acceleration.fs/2;
% %     temp.ttt=1:Settings.window*Settings.overlap_ratio:size(temp.ps,2);
% %% Calculate and plot FFT's of acceleration to determine proper higpass cutoff frequency
% % FFT pre-processing
% temp.window = Acceleration.fs * Settings.window;
% temp.noverlap = temp.window * Settings.overlap_ratio;
% temp.nfft = temp.window;
% 
% temp.points={Acceleration.awake_points,Acceleration.sleep_points};
% 
% % calculate psd with spectrogram with points
% figure('Name',[Mouse,' Normalized Spectrogram'])
% for j=1:2 % state
% for jj=1:Acceleration.nDim % Dimension
%     
%     temp.nPoints=size(temp.points{j},1);
%     temp.ps=cell(temp.nPoints,1);
%     for jjj=1:temp.nPoints
% 
%     [~,~,~,temp.ps{jjj}]=spectrogram(Acceleration.XYZ(temp.points{j}(jjj,1):temp.points{j}(jjj,2),jj),...
%     temp.window, temp.noverlap, temp.nfft, Acceleration.fs);
% 
%     end
%     temp.ps=cat(2,temp.ps{:});
%     
% % plot normalized (0-maxFreq Hz) spectrogram
% subplot(Acceleration.nDim,2,j+(jj-1)*2)
% temp.fff=0:Acceleration.fs/temp.nfft:Acceleration.fs/2;
% temp.ttt=1:Settings.window*Settings.overlap_ratio:size(temp.ps,2);
% temp.logical_fff=temp.fff<=Settings.maxFreq_acc;
% temp.ps_norm=(temp.ps(temp.logical_fff,:)./sum(temp.ps(temp.logical_fff,:)))*100; %normalized power spectrum with frequencies of interest (%)
% imagesc(temp.ttt,temp.fff(temp.logical_fff),temp.ps_norm)
% set(gca,'YDir','normal')
% colormap('PARULA')
% caxis([0,mean(max(temp.ps_norm))])
% colorbar              
% xlabel('Time (s)'),ylabel('Acceleration PSD (%)')
% 
% if j==1
%     title(['Awake CH',num2str(jj)])
% else
%     title(['Sleep CH',num2str(jj)]);
% end
%     
% Acceleration.PSD_cell{j,jj}=mean(temp.ps,2); % rows for awake/sleep, columns for channels
% 
% end
% end
% 
% clearvars temp
% 
% %% Some unnecessary plots and particularly options for awake/sleep distinguishment
% 
% Settings.LP_acc = .1; % Hz
%    
% % approximate dynamic acceleration per sampling time by differencing USE norm(X,Y,Z)
% Acceleration.dyn=(abs(diff(Acceleration.XYZ(:,1)))+...
%                   abs(diff(Acceleration.XYZ(:,2)))+... 
%                   abs(diff(Acceleration.XYZ(:,3))))*Acceleration.fs; % m/s^2[*(1/g) normalized?
%                             
% % approximate XYZ dynamic acceleration by differencing 
% Acceleration.XYZ_difference=[diff(Acceleration.XYZ(:,1)),...
%                              diff(Acceleration.XYZ(:,2)),...
%                              diff(Acceleration.XYZ(:,3))];
%                          figure('Name','Shows XYZ of differenced Acceleration')
%                          plot(Acceleration.t(1:end-1),Acceleration.XYZ_difference)
% 
% % approximate XYZ static acceleration by applying a low-pass filter
% Acceleration.XYZ_stat=[LowPassfilter(2,Settings.LP_acc,Acceleration.fs,Acceleration.XYZ(:,1)),...
%                        LowPassfilter(2,Settings.LP_acc,Acceleration.fs,Acceleration.XYZ(:,2)),...
%                        LowPassfilter(2,Settings.LP_acc,Acceleration.fs,Acceleration.XYZ(:,3))];
%                    figure('Name','Shows XYZ of filtered Acceleration (Low-pass 1Hz)')
%                    plot(Acceleration.t,Acceleration.XYZ_stat)
%                    
% figure('Name','Shows XYZ of filtered Acceleration (High-pass 1Hz)')
% plot(Acceleration.t,Acceleration.XYZ_dyn)
% 
% figure('Name','Shows length of XYZ dynamic acceleration vector (approximated by High-pass filter (1 Hz))')
% plot(Acceleration.t,Acceleration.dyn)
% 
% %% Plot spectrograms with standard settings with no distinguishment between sleep and awake
% % FFT pre-processing
% temp.window = Electrical.fs * Settings.window;
% temp.noverlap = temp.window * Settings.overlap_ratio;
% temp.nfft = temp.window;
% 
% figure('Name','Standard spectrograms for each channel')
% for j=1:Electrical.nChannel
%     subplot(Electrical.nChannel,1,j)
%     spectrogram(Electrical.CH1234(:,j),...
%     temp.window, temp.noverlap, temp.nfft, Electrical.fs,'yaxis');
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %%
% for i=1:N
%     clearvars -except i, close all
% N=1e3;
% x=randn(N,2);
% y=randn(2,N);
% 
% tic
% x=mean(x');
% time(i,1)=toc;
% 
% tic
% y=mean(y);
% time(i,2)=toc;
% end
% 
% %%
% for i=1:N
% N=1e4;
% 
% tic
% x=randn(1,N);
% A_time(i,2)=toc;
% 
% clear all
% 
% tic
% x=randn(N,1);
% A_time(i,1)=toc;
% 
% clear all
% 
% end
% %%
% 
% val=dyn_awake_threshold_new(2);
% pdf_Awake_rel=sym(@(d)1./(1+pdf_Sleep(d)./pdf_Awake(d)));
% pdf_Sleep_rel=sym(@(d)1./(1+pdf_Awake(d)./pdf_Sleep(d)));
% figure
% plot(dyn,pdf_Awake_rel(dyn),'r'),hold on
% plot(dyn,pdf_Sleep_rel(dyn),'b')
% plot(dyn,pdf_Awake(dyn),'-r')
% plot(dyn,pdf_Sleep(dyn),'-b')
% 
% 1/(1+pdf_Sleep(dyn_awake_threshold_new(2))/pdf_Awake(dyn_awake_threshold_new(2)))
% 
% 
% state012_select_sleep{j}=state012{j}(:,Acceleration.dyn_mean<dyn_sleep_threshold&Acceleration.dyn_mean<dyn_awake_threshold);
% state012{j}(:,Acceleration.dyn_mean<dyn_sleep_threshold&Acceleration.dyn_mean<dyn_awake_threshold)=ones(nPoints(j),size(state012_select_sleep{j},2))*index(2);
% 
% state012_select_awake{j}=state012{j}(:,Acceleration.dyn_mean>dyn_sleep_threshold&Acceleration.dyn_mean>dyn_awake_threshold);
% %% Acceleration state plot    
% figure('Name',[Mouse,' Acceleration state'])
% for j=1:size(points{1},1)
%     for jj=1:2:size(points{1}{j},2)
%     for jjj=1:size(Acceleration.XYZ,1)
%     
%     Acceleration.XYZ_transpose=Acceleration.XYZ';
%     plot(Acceleration.t(points{1}{j}(jj):points{1}{j}(jj+1)),...
%          Acceleration.XYZ_transpose(points{1}{j}(jj):points{1}{j}(jj+1),jjj),...
%          'Color',c(j,:),'LineStyle',ls{jjj}), hold on
% 
%     end
%     end
% end
% xlabel('Time (s)'), ylabel('Acceleration ^{1}/_{g}')
% %% old
% if ~strcmp(FileName{experiment_index},'10631-02_wt_D1.mat')     
% % truncate data 30 min each side
% % truncate acceleration data
% logical_trunc_acc=Acceleration.t>30*60&Acceleration.t<Acceleration.t(end)-1/Acceleration.fs-30*60;
% Acceleration.v=Acceleration.v(logical_trunc_acc);
% Acceleration.XYZ=Acceleration.XYZ(:,logical_trunc_acc);
% Acceleration.t=Acceleration.t(logical_trunc_acc);
% % truncate electrical data
% logical_trunc_ele=Electrical.t>30*60&Electrical.t<Electrical.t(end)-1/Electrical.fs-30*60;
% Electrical.CH1234=Electrical.CH1234(:,logical_trunc_ele);
% Electrical.t=Electrical.t(logical_trunc_ele);
% else
% end
%     
% logical_trunc_acc=Acceleration.t>60&Acceleration.t<Acceleration.t(end)-1/Acceleration.fs-60;
% logical_trunc_ele=Electrical.t>60&Electrical.t<Electrical.t(end)-1/Electrical.fs-60;
% %%
% % calculate absolute change of acceleration
% Position.v=sqrt(diff(Position.XYZ(1,:)).^2  ...
%                +diff(Position.XYZ(2,:)).^2+ ... 
%                 diff(Position.XYZ(3,:)).^2)*Position.fs;
% 
% %% plot 3d trace of mice path in cage for a check of coordinate interpretation -nope
% figure('Name','3d trace of mouse position')
% scatter3(Position.XYZ(1,logical120{1}(1,:)),Position.XYZ(2,logical120{1}(1,:)),Position.XYZ(3,logical120{1}(1,:))), hold on
% scatter3(Position.XYZ(1,logical120{1}(2,:)),Position.XYZ(2,logical120{1}(2,:)),Position.XYZ(3,logical120{1}(2,:)))
% %%
% % load('10425-01_het','Electrical','Position')  
% 
% 
% spectrogram(dataLP_select,...
%     paramsspect.window, paramsspect.noverlap, paramsspect.nfft, Settings.fs, 'yaxis');
% 
% plot(v,pdf(gmd,v'),'LineWidth',2.5)
% 
%    % Position.XYZ(:,7.1714e6+1)=[0;0;0] but trunc_position=7.1714e6;
%    
% %%%%%%%%%% this is all only relevant after last step of previous
% % kstest2 to see if statistically signigicant difference between
% % distributions
% 
% 
% % plot histogram to see if distribution (histcounts and hist)
% figure('Name','Histogram of speed bins')
% histogram(Position.v)
% 
% clearvars, close all
% tic
% load('10425-01_het','Electrical','Position')
% time(1)=toc;
% %%
% tic
% 
% Electrical.t=0:1/Electrical.fs:(length(Electrical.CH1234)-1)/Electrical.fs; % \muV
% Position.t=0:1/Position.fs:(length(Position.XYZ)-1)/Position.fs; % \ m? --> some kind of displacement or coordinates? 
%    
% Position.v=sqrt(diff(Position.XYZ(1,:)).^2+diff(Position.XYZ(2,:)).^2+diff(Position.XYZ(2,:)).^2)*Position.fs;
%    
% % truncate v and t 30 min each side
% Position.v=Position.v(Position.t>30*60&Position.t<Position.t(end)-1/Position.fs-30*60);
%    
% % additional small truncation with so that we have 1 min bins
% % timeBin=20; % s
% 
% timeBin=1:1:1000;
% timeBin=[.01,.02,.05,.1,.5,timeBin,1005:5:2000,2050:50:3600];
% for ii=1:length(timeBin)
% % ii=1;
% 
% Position.v_calc=Position.v(1:end-mod(length(Position.v),timeBin(ii)*Position.fs));
%    
% % reshape v into even parts (1 min bins)
% Position.v_calc=reshape(Position.v_calc,[timeBin(ii)*Position.fs,length(Position.v_calc)/(timeBin(ii)*Position.fs)]);
% 
% % mean reshaped v
% Position.v_calc=mean(Position.v_calc,1);
%    
% % fit mixture of gaussians with fitgmdist(data,2)
% gmd=fitgmdist(Position.v_calc',2);
% 
% % Indexing: index(1) is the index integer for the distribution with the higher mean
% % which should be awake, alseep is then index(2) with the lower mean
% [~,index]=sort(gmd.mu,'descend');
% 
% % compute individual and combined probability density function
% % Normalized from -inf to inf
% pd_Awake=@(v)1/sqrt(2*pi*gmd.Sigma(index(1)))*exp(-((v-gmd.mu(index(1))).^2)/(2*gmd.Sigma(index(1))));
% pd_Asleep=@(v)1/sqrt(2*pi*gmd.Sigma(index(2)))*exp(-((v-gmd.mu(index(2))).^2)/(2*gmd.Sigma(index(2))));
% PDF=@(v)gmd.ComponentProportion(index(1))*pd_Awake(v)+gmd.ComponentProportion(index(2))*pd_Asleep(v);
% % Normalized from 0 to inf
% % pd_Awake= @(v)sqrt(pi*gmd.Sigma(index(1))/2)*(1+erf(gmd.mu(index(1))/sqrt(2*gmd.Sigma(index(1)))))...
% %              *exp(-((v-gmd.mu(index(1))).^2)/(2*gmd.Sigma(index(1))));
% % pd_Asleep=@(v)sqrt(pi*gmd.Sigma(index(1))/2)*(1+erf(gmd.mu(index(1))/sqrt(2*gmd.Sigma(index(1)))))...
% %              *exp(-((v-gmd.mu(index(1))).^2)/(2*gmd.Sigma(index(1))));
% % PDF=@(v)gmd.ComponentProportion(index(1))*pd_Awake(v)+gmd.ComponentProportion(index(2))*pd_Asleep(v);
% 
% v=-1000:0.01:1000;
% % figure('Name','pd_Awake')
% % plot(v,pd_Awake(v)*0.001)
% % figure('Name','pd_Asleep')
% % plot(v,pd_Asleep(v)*0.001)
% % figure('Name','pd')
% % plot(v,pdf(gmd,v'))
% 
% PDF_VAL=PDF(v)*0.01;
% errorish(ii)=sum(PDF_VAL(v<0))*100;
% iterations(ii)=gmd.NumIterations;
% 
% % figure('Name','Histogram with pd')
% % histogram(Position.v,200), hold on
% % plot(v,pdf(gmd,v')*length(Position.v),'Marker','x')
% % plot(v,PDF(v)*length(Position.v),'Marker','o')
% % xlim([0,100])
%    
% 
% 
% % figure out what these are as probability densities and decide for which speed THRESHOLD bin is likely to belong to where 
%    
% time(2,ii)=toc;
% 
% end
% 
% figure
% plot(timeBin(1:ii-1),errorish), hold on
% plot(timeBin(1:ii-1),iterations)
% xlabel('time length of time bin'), ylabel('%error and #iterations')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sleep compiler 


clearvars
close all

% EF = uigetdir('','Select Experiment Folder (EF)');
% cd(EF)
% E=dir('*_Analyzed.mat');
% nE=length(E);

% Select to be analyzed .mat files
[E_Name,PathName] = uigetfile('*Analyzed.mat','Select the file to analyse','MultiSelect', 'on');
if ischar(E_Name) %if-statement for the case that only one file is selected
    E_Name={E_Name}; 
end
cd(PathName)
nE=length(E_Name);
for i=1:nE
    E(i).name=E_Name{i};
end

load('DataTable_StartTime')

% Interpolation: When the state is unknown interpolate to surrounding state
% if the length of unknown segment is also under appropriate defined interpolation time.
AwakeInterpolationTime=3*60;
 NREMInterpolationTime=3*60;
  REMInterpolationTime=30;

% Rule 1: QA between REM or NREM and REM is likely also REM (<REMInterpolationTime)
REM2QAInterpolationTime=30;
% Rule 2: REM does not follow Awake(=AA||QA), instead this is probably QA. Only AA should be used (QA untrustworthy). 
REMafterAwakeCutoffTime=3*60;

% Get sampling rate and assume to it be constant across measurements
load(E(1).name,'Electrical')
fs=Electrical.fs;
clearvars Electrical

c=colormap(lines);
c([2,1],:)=c([1,2],:);
close(gcf)

% Segment length for circadian statistics
segmentLength=round(60*60*fs);
segmentOverlap=round(segmentLength*(9/10));
dt=(segmentLength-segmentOverlap)/fs;

% Define between which hours it is light and dark
Light=[7,19];
% Dark=[19,7];
% maxTime=96; % hours
stateNames={'Awake';'NREM';'REM'};
nstateNames=length(stateNames);
% cut_off=0.5*10^-4;
cut_off=2.5*10^-5;

% nE=2;
circadianTime=cell(nE,1);
circadianStatistics=cell(nE,1);
StatePointsCell=cell(nE,1);
circadianFiltCell=cell(nE,1);
circadianFiltTime=cell(nE,1);
fs_circadian=zeros(nE,1);
for i=1:nE
    
    load(E(i).name,'Electrical_P')
    nanLogical=~isnan(Electrical_P(:,1));
    [~,state_index]=max(Electrical_P(nanLogical,:),[],2);
    clearvars Electrical_P
    nanLogical_double=double(nanLogical);
    nanLogical_double(nanLogical)=state_index;
    state_index=nanLogical_double;
    clearvars nanLogical_double nanLogical
    
    nElectrical=size(state_index,1);
    
    Unknown=logical2points(state_index==0);
    AA=logical2points(state_index==1);
    NREM=logical2points(state_index==2);
    REM=logical2points(state_index==3);
    QA=logical2points(state_index==4);
    clearvars state_index
    
    nUnknown=size(Unknown,1); % Interpolation
    nQA=size(QA,1); % Rule 1
    nREM=size(REM,1); % Rule 2
    nAA=size(AA,1);
    nNREM=size(NREM,1);
    
    StateID=[repmat({'Unknown'},[nUnknown,1]);...
             repmat({'AA'},[nAA,1]);...
             repmat({'NREM'},[nNREM,1]);...
             repmat({'REM'},[nREM,1]);...
             repmat({'QA'},[nQA,1])];
         
    StatePoints=table([Unknown;AA;NREM;REM;QA],StateID,'variablenames',{'Points','State'});
    StatePoints=sortrows(StatePoints,'Points','ascend');
    
%     % test table, interpolation and rule
%     StatePoints=table([1,100;101,200;201,300;301,400;401,500;501,600;601,700;701,800;801,900;901,1000;1001,1100;1101,1200;1201,1300;1301,1400;1401,1500;1501,1600;1601,1700;1701,1800;1801,1900;1901,2000;2001,2100;2101,2200;2201,2300;2301,2400;2401,2500;2501,2600;2601,2700;2701,2800;2801,2900;2901,3000;3001,3100;3101,3200;3201,3300],...
%     {'Unknown','AA','Unknown','AA','QA','Unknown','QA','NREM','Unknown','NREM','REM','Unknown','REM','QA','REM','AA','REM','AA','Unknown','NREM','Unknown','QA','REM','Unknown','AA','Unknown','REM','NREM','Unknown','REM','Unknown','NREM','Unknown'}','variablenames',{'Points','State'});
%     StatePoints=table([1,1;2,2;3,round(AwakeInterpolationTime*fs)+3;round(AwakeInterpolationTime*fs)+4,round(AwakeInterpolationTime*fs)+4;round(AwakeInterpolationTime*fs)+5,round(AwakeInterpolationTime*fs)+5],...
%                 {'Unknown','AA','Unknown',                              'AA',                                                                   'Unknown'}','variablenames',{'Points','State'});
%     
%     nUnknown=sum(strcmp(StatePoints.State,'Unknown')); % Interpolation
%     nQA=sum(strcmp(StatePoints.State,'QA')); % Rule 1
%     nREM=sum(strcmp(StatePoints.State,'REM')); % Rule 2
    
    nSegment=size(StatePoints.Points,1);
            
    % Interpolation, maybe split up points in halfs and assign state
    UnknownIndex=find(strcmp(StatePoints.State,'Unknown'));
    nInterp=0;
    for ii=1:nUnknown
        if UnknownIndex(ii)~=1&&UnknownIndex(ii)~=nSegment
            % AA/QA-AA/QA (1,1)
            unknownLength=(StatePoints.Points(UnknownIndex(ii),2)-StatePoints.Points(UnknownIndex(ii),1)+1);
            if (strcmp(StatePoints.State{UnknownIndex(ii)-1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)-1},'QA'))&&...
               (strcmp(StatePoints.State{UnknownIndex(ii)+1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)+1},'QA'))
                if unknownLength/fs<AwakeInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'AA'}; % ignoring QA for now
                    nInterp=nInterp+1;
                end
            % NREM-NREM (2,2)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'NREM')&&...
                   strcmp(StatePoints.State{UnknownIndex(ii)+1},'NREM')
                if unknownLength/fs<NREMInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'NREM'};
                    nInterp=nInterp+1;
                end
            % REM-REM (3,3)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'REM')&&...
                   strcmp(StatePoints.State{UnknownIndex(ii)+1},'REM')
                if unknownLength/fs<REMInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'REM'};
                    nInterp=nInterp+1;
                end
            % AA/QA-NREM (1,2)
            elseif (strcmp(StatePoints.State{UnknownIndex(ii)-1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)-1},'QA'))&&...
                    strcmp(StatePoints.State{UnknownIndex(ii)+1},'NREM')
                if unknownLength/fs<(AwakeInterpolationTime+NREMInterpolationTime)/4
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp=nInterp+1;
                end
            % NREM-AA/QA (2,1)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'NREM')&&...
                  (strcmp(StatePoints.State{UnknownIndex(ii)+1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)+1},'QA'))
                if unknownLength/fs<(AwakeInterpolationTime+NREMInterpolationTime)/4
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp=nInterp+1;
                end
            % AA/QA-REM (1,3)
            elseif (strcmp(StatePoints.State{UnknownIndex(ii)-1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)-1},'QA'))&&...
                   strcmp(StatePoints.State{UnknownIndex(ii)+1},'REM')
                if unknownLength/fs<REMInterpolationTime/2
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp=nInterp+1;
                elseif unknownLength/fs<AwakeInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'AA'}; % ignoring QA
                    nInterp=nInterp+1;
                end
            % REM-AA/QA (3,1)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'REM')&&...
                   (strcmp(StatePoints.State{UnknownIndex(ii)+1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)+1},'QA'))
                if unknownLength/fs<REMInterpolationTime/2
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp=nInterp+1;
                elseif unknownLength/fs<AwakeInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'AA'}; % ignoring QA
                    nInterp=nInterp+1;
                end
            % NREM-REM (2,3)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'NREM')&&...
                   strcmp(StatePoints.State{UnknownIndex(ii)+1},'REM')
                if unknownLength/fs<REMInterpolationTime/2
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp=nInterp+1;
                elseif unknownLength/fs<NREMInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'NREM'}; % ignoring QA
                    nInterp=nInterp+1;
                end
            % NREM-REM (3,2)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'REM')&&...
                   strcmp(StatePoints.State{UnknownIndex(ii)+1},'NREM')
                if unknownLength/fs<REMInterpolationTime/2
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp=nInterp+1;
                elseif unknownLength/fs<NREMInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'NREM'}; % ignoring QA
                    nInterp=nInterp+1;
                end
            end
        end
    end
    StatePoints(strcmp(StatePoints.State,'Remove'),:)=[];
    nSegment=size(StatePoints.Points,1);
    nAA=sum(strcmp(StatePoints.State,'AA'));
    nNREM=sum(strcmp(StatePoints.State,'NREM'));
    nREM=sum(strcmp(StatePoints.State,'REM'));
    nQA=sum(strcmp(StatePoints.State,'QA'));
    nUnknown=sum(strcmp(StatePoints.State,'Unknown'));
    
    % Rule 1
    QAIndex=find(strcmp(StatePoints.State,'QA'));
    nQAInterp=0;
    for ii=1:nQA
        if (strcmp(StatePoints.State{QAIndex(ii)-1},'REM')||strcmp(StatePoints.State{QAIndex(ii)-1},'NREM'))&&...
            strcmp(StatePoints.State{QAIndex(ii)+1},'REM')
            if (StatePoints.Points(QAIndex(ii),2)-StatePoints.Points(QAIndex(ii),1)+1)/fs<REM2QAInterpolationTime
                StatePoints.State(QAIndex(ii))={'REM'};
                nQAInterp=nQAInterp+1;
            end
        end
    end
    
    % Rule 2
    REMIndex=find(strcmp(StatePoints.State,'REM'));
    nREM2QA=0;
%     figure
    for ii=1:nREM
        if REMIndex(ii)~=1
            if strcmp(StatePoints.State{REMIndex(ii)-1},'AA')
                diffREM=(StatePoints.Points(REMIndex(ii),1)-StatePoints.Points(strcmp(StatePoints.State,'AA'),2)+1)/fs;
                if any((diffREM<REMafterAwakeCutoffTime)&(diffREM>0))
                    StatePoints.State(REMIndex(ii))={'QA'};
                    nREM2QA=nREM2QA+1;
                end
%                 plot(diffREM),hold on
            end
        end
    end
    
    % Merge AA and QA to Awake
    StatePoints.State(strcmp(StatePoints.State,'AA')|strcmp(StatePoints.State,'QA'))={'Awake'};
    % Calculate statistics (maybe exclude segments bordered by Unknown)
    
    % Anova
    % Plot boxplots
    % pre
    UnknownLog=points2logical(StatePoints.Points(strcmp(StatePoints.State,'Unknown'),:),nElectrical);
    AwakeLog  =points2logical(StatePoints.Points(strcmp(StatePoints.State,'Awake'),:),nElectrical);
    NREMLog   =points2logical(StatePoints.Points(strcmp(StatePoints.State,'NREM'),:) ,nElectrical);
    REMLog    =points2logical(StatePoints.Points(strcmp(StatePoints.State,'REM'),:)  ,nElectrical);
    UnknownPoints=logical2points(UnknownLog); clearvars UnknownLog
    AwakePoints  =logical2points(AwakeLog);
    NREMPoints   =logical2points(NREMLog);
    REMPoints    =logical2points(REMLog);
    
    % Adjust and store StatePoints table
    StateID=[repmat({'Unknown'},[size(UnknownPoints,1),1]);...
             repmat({'Awake'},[size(AwakePoints,1),1]);...
             repmat({'NREM'},[size(NREMPoints,1),1]);...
             repmat({'REM'},[size(REMPoints,1),1])];
    StatePoints=table([UnknownPoints;AwakePoints;NREMPoints;REMPoints],StateID,'variablenames',{'Points','State'});
    StatePoints=sortrows(StatePoints,'Points','ascend');
    StatePointsCell{i}=StatePoints;
    
%     %
%     startTime=DataTable_StartTime.StartTime(strcmp(DataTable_StartTime.Name,E(i).name));
%     startTime=hour(startTime)*60*60+minute(startTime)*60+second(startTime);
%     t=mod(Electrical.t+startTime,24*60*60);
%     Light_logical=t>Light(1)*60*60&t<Light(2)*60*60;
%     Dark_logical=~Light_Logical;
%     Light_logical=Light_logical&points2logical(Electrical.points,Electrical.n);
%     Dark_logical=Dark_logical&points2logical(Electrical.points,Electrical.n);
    
    % Calculate circadian statistics
    KnownPoints=logical2points(~points2logical(StatePoints.Points(strcmp(StatePoints.State,'Unknown'),:),nElectrical));   
    [AwakeLog,midPoint]=buffer_dis(AwakeLog,segmentLength,segmentOverlap,KnownPoints); AwakeLog=sum(AwakeLog,1)';
    [NREMLog,~]        =buffer_dis(NREMLog ,segmentLength,segmentOverlap,KnownPoints); NREMLog=sum(NREMLog,1)';
    [REMLog,~]         =buffer_dis(REMLog  ,segmentLength,segmentOverlap,KnownPoints); REMLog=sum(REMLog,1)';
    
    StartTime=DataTable_StartTime.StartTime(strcmp(DataTable_StartTime.Name,E(i).name));
    circadianTime{i}=midPoint/fs+hour(StartTime)*3600+minute(StartTime)*60+second(StartTime);
    nmidPoint=length(midPoint);
    circadianStatistics{i}=zeros(nmidPoint,3);
    circadianStatistics{i}=[AwakeLog,NREMLog,REMLog];
    circadianStatistics{i}(:,1:3)=(circadianStatistics{i}(:,1:3)./repmat(sum(circadianStatistics{i}(:,1:3),2),[1,3]))*100;
    circadianStatistics{i}(:,3)=(circadianStatistics{i}(:,3)./sum(circadianStatistics{i}(:,2:3),2))*100;
    circadianStatistics{i}(isnan(circadianStatistics{i}(:,3)),3)=0;
    
    % Filter
    fs_circadian=1/dt;
    circadianFiltTime{i}=(circadianTime{i}(1):dt:circadianTime{i}(end))';
    circadianFilt=zeros(length(circadianFiltTime{i}),3);
    for ii=1:nstateNames
        %%%%%%%%%%%%%%%%% Would also interpolate >24 hour periods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        circadianFilt(:,ii)=interp1(circadianTime{i},circadianStatistics{i}(:,ii),circadianFiltTime{i});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        circadianFilt(:,ii)=LowPassfilter(2, cut_off, fs_circadian, circadianFilt(:,ii));
    end
    circadianFiltCell{i}=circadianFilt;
    
end

%%
UpDi=triu(true(nE),1);
LoDi=tril(true(nE),-1);
nMaxlag=round(12*60*60*fs_circadian);
lagMatrix=zeros(nE,nE,nstateNames);
% lagMatrix2=zeros(nE,nE,nstateNames);
% lagMatrix3=zeros(nE,nE,nstateNames);
startTimeShift=zeros(nE,nE,nstateNames);
lagVec=zeros(nE,nE);
% lagVec2=zeros(nE,nstateNames);
CCMatrix=zeros(nE,nE,nstateNames);

% locationLocMin=zeros(nE,nstateNames);
% locationLocMax=zeros(nE,nstateNames);

% lagMatrix_save=cell(nstateNames,1);
% lagVecSave=cell(nstateNames,1);
% quitIterTime=1/10; % hours
for iii=1:nstateNames
    
%     for i=1:nE
%         locationLocMin(i,iii)=min(circadianFiltTime{i}(islocalmin(circadianFiltCell{i}(:,iii))))/3600;
%         locationLocMax(i,iii)=min(circadianFiltTime{i}(islocalmax(circadianFiltCell{i}(:,iii))))/3600;
%     end
    
    % calculate cross correlation
    for i=1:nE
    for ii=1:nE
    if UpDi(i,ii)
    % Find time delay between mice circadian rhythms
%     nWindow=length(circadianFiltTime{i});
%     circadianTempTime=(circadianFiltTime{i}(1)-3600*12:1/fs_circadian:circadianFiltTime{i}(end)+3600*12)';
%     circadianTempFilt=[zeros(floor((length(circadianTempTime)-length(circadianFiltCell{i}))/2),1);circadianFiltCell{i}(:,iii);zeros(ceil((length(circadianTempTime)-length(circadianFiltCell{i}))/2),1)];
%     referenceWave=sin((circadianTempTime/3600-((circadianTempTime(i))/3600+2))*(pi/12));
%     if iii==2
%         referenceWave=-referenceWave;
%     end
    nWindow=max([length(circadianFiltTime{i}),length(circadianFiltTime{ii})]);
    [C,~,L] = corrgram2(circadianFiltCell{i}(:,iii)-mean(circadianFiltCell{i}(:,iii)),circadianFiltCell{ii}(:,iii)-mean(circadianFiltCell{ii}(:,iii)),nMaxlag,nWindow,0,fs_circadian);

%     diffLoc=[locationLocMin(ii,iii)-locationLocMin(i,iii),locationLocMax(ii,iii)-locationLocMax(i,iii)];
%     [~,index]=min(abs(diffLoc));
%     lagMatrix3(i,ii,iii)=diffLoc(index);
    
%     if max(C(:))>1
%         break
%     end

    L=L/3600; % hours
    startTimeShift(i,ii,iii)=(circadianFiltTime{ii}(1)-circadianFiltTime{i}(1))/3600; % start time shift in hours
    [pks,locs,~,~]=findpeaks(C,L,'WidthReference','halfheight');
    if ~isempty(pks)
        [~,peak_closestto0_index]=min(abs(locs));
        lagMatrix(i,ii,iii)= locs(peak_closestto0_index); % shift i by this value to be aligned to ii
        CCMatrix(i,ii,iii) = pks (peak_closestto0_index);
    else
        lagMatrix(i,ii,iii)= L(nMaxlag+1);
        CCMatrix(i,ii,iii) = C(nMaxlag+1);
    end
%     lagMatrix2(i,ii,iii)=finddelay(circadianFiltCell{i}(:,iii),circadianFiltCell{ii}(:,iii),nMaxlag*2)/(fs_circadian*3600);
%     if abs(round(lagMatrix2(i,ii,iii)))>=12
%         lagMatrix2(i,ii,iii)=0;
%     end
    lagMatrix(i,ii,iii) =lagMatrix(i,ii,iii) +startTimeShift(i,ii,iii);
%     lagMatrix2(i,ii,iii)=lagMatrix2(i,ii,iii)+startTimeShift(i,ii,iii);
    
    
%         %%%%%%%%%%%%%%%%%%%%%%%
%     figure('name',[num2str(lagMatrix(i,ii,iii)),',',stateNames{iii},',',E(i).name])
%     subplot(1,2,1)
%     plot(circadianTempTime/3600,circadianTempFilt/max(circadianTempFilt)), hold on
%     plot(circadianTempTime/3600,referenceWave/max(referenceWave))
%     xtickVal=0:6:4*24;
% xticks(xtickVal)
% xticklabels(mod(xtickVal,24))
% xlim([min(startPlot),max(endPlot)])
% XLIM=get(gca,'XLim');
% YLIM=get(gca,'YLim');
% N=1e3;
% contourT=linspace(-11,109,1e3);
% contourP=YLIM(2)*ones(1,N);
% contourZ=zeros(1,N);
% colormap(gray)
% contourC=repmat([linspace(1,0,N/10),linspace(0,1,N/10)],[1,5]);
% surface([contourT;contourT],[contourP;contourP],[contourZ;contourZ],[contourC;contourC],...
%         'facecol','no',...
%         'edgecol','interp',...
%         'linew',5);
%     subplot(1,2,2)
%     plot(L,C);
%     %%%%%%%%%%%%%%%%%%%%%%%
    
    
    end    
    end
    end
    

    

    
%     [~,index]=min(sum(abs(lagMatrix(:,:,iii)),2)/(nE-1));
%     lagVec(:,iii)=lagMatrix(:,index,iii);
%     [~,index2]=min(sum(abs(lagMatrix2(:,:,iii)),2)/(nE-1));
%     lagVec2(:,iii)=lagMatrix2(:,index2,iii);
    
%     [lagVec2(:,1),lagMatrix2(:,:,1)]
%     nIter=0; index2=0;
%     while max(abs(lagMatrix{iii}(:)))>quitIterTime
%         lagVecAbs=sum(abs(lagMatrix{iii}),2)/(nE-1);
%         [~,index]=max(lagVecAbs);
%         if index==index2
%             index=randi(nE);
%         end
%         lagVec=max(abs(lagMatrix{iii}),[],2);
%         lagVal=lagVec(index);
%         lagMatrix{iii}(index,:)=lagMatrix{iii}(index,:)-lagVal;
%         lagMatrix{iii}(:,index)=lagMatrix{iii}(:,index)+lagVal;
%         lagMatrix{iii}(index,index)=0;
%         lagVecSave{iii}(index,1)=lagVecSave{iii}(index,1)-lagVal;
%         nIter=nIter+1; index2=index;
%     end
    
end
% pick row with lowest absolute delay on average
% align all other waves to this wave (column contains necessary shifts)

for i=1:nE
for ii=1:nE
if UpDi(i,ii)
    if CCMatrix(i,ii,1)>=CCMatrix(i,ii,2)
        lagVec(i,ii)=lagMatrix(i,ii,1);
    else
        lagVec(i,ii)=lagMatrix(i,ii,2);
    end
end
end
end
for i=1:nE
for ii=1:nE
    if LoDi(i,ii)
        lagVec(i,ii)=-lagVec(ii,i);
    end
end
end
[~,index]=min(sum(abs(lagVec),2)/(nE-1));
lagVec=lagVec(index,:)'; % adjust starting times with this vector (hours)

% lagVecSave=(lagVecSave{1}+lagVecSave{2})/2;
% lagVec11=( lagVec(:,1)+ lagVec(:,2))/2;
% lagVec22=(lagVec2(:,1)+lagVec2(:,2))/2;
% for i=1:nE
%     circadianFiltTime{i}=circadianFiltTime{i}-3600*lagVec(i);
% end
% lagMatrix
% lagMatrix2

%% use adjusted starting times and exclusion of border segments to calculate new boxplot statistics

%% Plot circadian
for i=1:nE
    
    figure
    for ii=1:nstateNames
        plot(circadianTime{i}/3600,circadianStatistics{i}(:,ii),'Color',c(ii,:)),hold on
    end
    legend(stateNames)
    xtickVal=0:6:4*24;
    xticks(xtickVal)
    xticklabels(mod(xtickVal,24))
    xlim([circadianTime{i}(1),circadianTime{i}(end)])

% figure
% dt=circadianTime{i}(2)-circadianTime{i}(1);
dt=mode(diff(circadianTime{i}));
% freq=-1/dt/2:1/(circadianTime{i}(end)-circadianTime{i}(1)):1/dt/2;
% if length(freq)~=length(circadianTime{i})
%     freq(end)=[];
% end
%     for ii=1:nstateNames
%         plot(freq,abs(fftshift(fft(circadianStatistics{i}(:,ii)))),'Color',c(ii,:)),hold on
%     end
% legend(stateNames)

figure
for ii=1:nstateNames
    plot(circadianFiltTime{i}/3600,circadianFiltCell{i}(:,ii),'Color',c(ii,:)),hold on
end
legend(stateNames)
xtickVal=0:6:4*24;
xticks(xtickVal)
xticklabels(mod(xtickVal,24))
xlim([circadianFiltTime{i}(1),circadianFiltTime{i}(end)]/3600)
ylim([0,100])

% figure
% dt=circadianTime{i}(2)-circadianTime{i}(1);
% freq=-1/dt/2:1/(circadianTime{i}(end)-circadianTime{i}(1)):1/dt/2;
% if length(freq)~=length(circadianTime{i})
%     freq(end)=[];
% end
%     for ii=1:nstateNames
%         plot(freq,abs(fftshift(fft(circadianFilt(:,ii)))),'Color',c(ii,:)),hold on
%     end
% legend(stateNames)

end

figure
for iii=1:nstateNames
subplot(1,3,iii)
startPlot=zeros(nE,1);
endPlot=zeros(nE,1);
for i=1:nE
    set(gcf,'name',num2str(cut_off))
    plot(circadianFiltTime{i}/3600-lagVec(i),circadianFiltCell{i}(:,iii)), hold on;
    startPlot(i)=circadianFiltTime{i}(1)/3600  -lagVec(i);
    endPlot(i)=  circadianFiltTime{i}(end)/3600-lagVec(i);
end
xtickVal=0:6:4*24;
xticks(xtickVal)
xticklabels(mod(xtickVal,24))
xlim([min(startPlot),max(endPlot)])
XLIM=get(gca,'XLim');
YLIM=get(gca,'YLim');
N=1e3;
contourT=linspace(-11,109,1e3);
contourP=YLIM(2)*ones(1,N);
contourZ=zeros(1,N);
colormap(gray)
contourC=repmat([linspace(1,0,N/10),linspace(0,1,N/10)],[1,5]);
surface([contourT;contourT],[contourP;contourP],[contourZ;contourZ],[contourC;contourC],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',5);
end
    
    
    
    
% ylim([0,100])

% 
% state_index2=state_index;
% state_index2(points2logical(StatePoints.Points(strcmp(StatePoints.State,'QA'),:),nElectrical))=4;
% state_index2(points2logical(StatePoints.Points(strcmp(StatePoints.State,'REM'),:),nElectrical))=3;
% state_index2(points2logical(StatePoints.Points(strcmp(StatePoints.State,'NREM'),:),nElectrical))=2;
% state_index2(points2logical(StatePoints.Points(strcmp(StatePoints.State,'AA'),:),nElectrical))=1;
% state_index2(points2logical(StatePoints.Points(strcmp(StatePoints.State,'Unknown'),:),nElectrical))=0;
