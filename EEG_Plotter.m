% Plot mean PSD for each mouse and genotype
clearvars, close all
%% Import
[FileName,PathName] = uigetfile('DataTable*.mat','Select the file to analyse','MultiSelect', 'on');
if ischar(FileName) %if-statement for the case that only one file is selected
    FileName={FileName}; 
end
cd(PathName)
nTables=length(FileName);

DataTable_cell=cell(nTables,1);
for i=1:nTables
    temp=load(FileName{i});
    DataTable_cell{i}=temp.DataTable;
end
DataTable=cat(1,DataTable_cell{:});

load('Settings')

%% Constants defined here (can be changed to input prompt)
maxFreq=10;
SEM_transparency=0.2;
normFreq=50;

% define fff by assuming fs to be constant
fs=DataTable.fs(1);
fff=0:fs/round(fs*Settings.window_FFT):fs/2;

% define channel names
Channel_Names(1,1)="right-S"; % CH1=right somatosensory cortex
Channel_Names(2,1)="right-M"; % CH2=right motor cortex
Channel_Names(3,1)="left-M";  % CH3= left motor cortex
Channel_Names(4,1)="left-S";  % CH4= left somatosensory cortex

% change channel names in table
[uniq.Channel,~,uniq.Channel_index]=unique(DataTable.Channel);
nChannel=length(uniq.Channel);
for i=1:nChannel
    DataTable.Channel(uniq.Channel_index==i)=Channel_Names(i,1);
end

% define frequency bands
delta=[1,4];
theta=[5,10];
beta=[13,30];
gamma=[30,50];
bandstr=["delta";"theta";"beta";"gamma"];
nBand=length(bandstr);

%% Plot settings
% plot property 'Color'
c=colormap(lines); % colors to most contrasting color map
% switch around colors (for wt and het)
c([1,2],:)=c([2,1],:);%%%%%%%%%%%%%%%
close(gcf)

% plot property 'LineStyle'
ls{1}='-';
ls{2}='--';
ls{3}=':';
ls{4}='-.';

% plot property 'LineWidth'
lw=1.5;

%% 

nState_full=length(unique(DataTable.State));

uniq.Name_full=unique(DataTable.Name);
uniq.State_full=unique(DataTable.State);

count=0;
uniq.Channel=unique(DataTable.Channel);
nChannel=length(uniq.Channel);
for i=1:nChannel
    logical.channel=strcmp(DataTable.Channel,uniq.Channel(i));
    
    figure('Name',[uniq.Channel{i},' Mean over Days'])
    legend_cell=cell(nState_full,2);
    handle_cell=cell(nState_full,2);
    
    uniq.Name=unique(DataTable.Name(logical.channel));
    
    for ii=1:length(uniq.Name)
    logical.channel_name=logical.channel&strcmp(DataTable.Name,uniq.Name{ii});
    uniq.State=unique(DataTable.State(logical.channel_name));
    
    for iii=1:length(uniq.State)
    logical.channel_name_state=logical.channel_name&strcmp(DataTable.State,uniq.State{iii});
        
    if sum(logical.channel_name_state)
    
    count=count+1;    
        
    % do mean over days and save in cell for temp table
    Power_Matrix=[DataTable.Power{logical.channel_name_state}];
    [temptable.Power{count,1},SEM_meanday]=meanSEM(Power_Matrix);
    
    % additional strings for temp table
    temptable.Channel(count,1)=uniq.Channel(i);
    temptable.Genotype(count,1)=DataTable.Genotype(find(logical.channel_name_state,1,'first'));
    temptable.State(count,1)=uniq.State(iii);
    temptable.Name(count,1)=uniq.Name(ii);
    TEMP(count)=sum(logical.channel_name_state);
    
    % name index code (line style consistent)
    name_index=find(strcmp(uniq.Name_full,uniq.Name(ii)));
    name_index_c=mod(name_index-1,size(c,1))+1;
    name_index_ls=mod(name_index-1,length(ls))+1;
    
    % state
    state_index=find(strcmp(uniq.State_full,uniq.State(iii)));
    
    % absolute plot
    subplot(2,2,(state_index-1)*2+1) % index for plot (2x2 : state x norm)
    h=plot(fff,temptable.Power{count},...
        'Color',c(name_index_c,:),'LineStyle',ls{name_index_ls},'LineWidth',lw); hold on
    patch([fff,...
           fff(end:-1:1),...
           fff(1)],...
          [temptable.Power{count}-SEM_meanday;...
           temptable.Power{count}(end:-1:1)+SEM_meanday(end:-1:1);...
           temptable.Power{count}(1)],...
           c(name_index_c,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
    title([uniq.Channel{i},' ',uniq.State_full{state_index},' -absolute'])
    xlim([Settings.HP_ele,maxFreq])
    legend_cell{iii,1}=[legend_cell{iii,1},DataTable.Name(find(logical.channel_name_state,1,'first'))];
    handle_cell{iii,1}=[handle_cell{iii,1},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    
    % normalization (%)
    Power_Matrix=(Power_Matrix./bandpower2(Power_Matrix,fff,[Settings.HP_ele,normFreq]))*100; 
    [Power_meanday_normalized,SEM_meanday_normalized]=meanSEM(Power_Matrix);
    
    % normalized plot
    subplot(2,2,state_index*2) % index for plot
    h=plot(fff,Power_meanday_normalized,...
        'Color',c(name_index_c,:),'LineStyle',ls{name_index_ls},'LineWidth',lw); hold on
    patch([fff,...
           fff(end:-1:1),...
           fff(1)],...
          [Power_meanday_normalized-SEM_meanday_normalized;...
           Power_meanday_normalized(end:-1:1)+SEM_meanday_normalized(end:-1:1);...
           Power_meanday_normalized(1)],...
           c(name_index_c,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
    title([uniq.Channel{i},' ',uniq.State_full{state_index},' -normalized'])
    xlim([Settings.HP_ele,maxFreq])
    legend_cell{iii,2}=[legend_cell{iii,2},DataTable.Name(find(logical.channel_name_state,1,'first'))];
    handle_cell{iii,2}=[handle_cell{iii,2},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (%)')
    
    end % end of if-statement
    
    end
    end

for l=1:nState_full
    for ll=1:2
        subplot(nState_full,2,ll+(l-1)*nState_full)
        legend(handle_cell{l,ll},legend_cell{l,ll})
    end
end

end 

% clearvars handle_cell legend_cell
%% Genotype over Names
% in and after this section fs is assumed to be constant

% pre
uniq.State_full=unique(temptable.State);
nState_full=length(uniq.State_full);
uniq.Genotype_full=unique(temptable.Genotype);
nGenotype_full=length(uniq.Genotype_full);

uniq.Channel=unique(temptable.Channel);
nChannel=length(uniq.Channel);
for i=1:nChannel
    logical.channel=strcmp(temptable.Channel,uniq.Channel(i));
    
    figure('Name',[uniq.Channel{i},' Mean over Genotype'])
    legend_cell=cell(nState_full,2);
    handle_cell=cell(nState_full,2);
    
    uniq.State=unique(temptable.State(logical.channel));
    nState=length(uniq.State);
    for ii=1:nState
    logical.channel_state=logical.channel&strcmp(temptable.State,uniq.State(ii));
    
    uniq.Genotype=unique(temptable.Genotype);
    nGenotype=length(uniq.Genotype);
    for iii=1:nGenotype
    logical.channel_state_genotype=logical.channel_state&strcmp(temptable.Genotype,uniq.Genotype(iii));
    
    % state index
    state_index=find(strcmp(uniq.State_full,uniq.State(ii)));
    % genotype index
    genotype_index=find(strcmp(uniq.Genotype_full,uniq.Genotype(iii)));
    
    % absolute
    Power_Matrix=cat(2,temptable.Power{logical.channel_state_genotype});
    [temp.power,temp.SEM]=meanSEM(Power_Matrix);
    
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
    title([uniq.Channel{i},' ',uniq.State_full{state_index},' -absolute'])
    xlim([Settings.HP_ele,maxFreq])
    legend_cell{ii,1}=[legend_cell{ii,1},temptable.Genotype(find(logical.channel_state_genotype,1,'first'))];
    handle_cell{ii,1}=[handle_cell{ii,1},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    
    % normalization (%)
    Power_Matrix=(Power_Matrix./bandpower2(Power_Matrix,fff,[Settings.HP_ele,normFreq]))*100; 
    [temp.power_norm,temp.SEM_norm]=meanSEM(Power_Matrix);
    
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
    title([uniq.Channel{i},' ',uniq.State_full{state_index},' -normalized'])
    xlim([Settings.HP_ele,maxFreq])
    legend_cell{ii,2}=[legend_cell{ii,2},temptable.Genotype(find(logical.channel_state_genotype,1,'first'))];
    handle_cell{ii,2}=[handle_cell{ii,2},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (%)')
    
    end
    end
    

for l=1:nState_full
    for ll=1:2
        subplot(nState_full,2,ll+(l-1)*nState_full)
        legend(handle_cell{l,ll},legend_cell{l,ll})
    end
end
    
end

% Average over areas
count=0;
temptable.Area=strings(length(temptable.Channel),1);
[uniq.Channel,~,uniq.Channel_index]=unique(temptable.Channel);
nChannel=length(uniq.Channel);
for i=1:nChannel
    temp.string=strsplit(uniq.Channel(i),'-');
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
    count=count+1;
    temptable2.Area(count,1)=temptable.Area(find(logical.area_state_name,1,'first'));
    temptable2.State(count,1)=temptable.State(find(logical.area_state_name,1,'first'));
    temptable2.Genotype(count,1)=temptable.Genotype(find(logical.area_state_name,1,'first'));
    temptable2.Power{count,1}=mean(cat(2,temptable.Power{logical.area_state_name}),2);
end

end
end
end       

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
bandtable.Power{count,1}=[bandpower2(temp.power_arr,fff,delta);...
                          bandpower2(temp.power_arr,fff,theta);...
                          bandpower2(temp.power_arr,fff,beta);...
                          bandpower2(temp.power_arr,fff,gamma);...
                          bandpower2(temp.power_arr,fff)];

        end
    end
    
    for l=1:nState_full
    for ll=1:2
        subplot(nState_full,2,ll+(l-1)*nState_full)
        legend(handle_cell{l,ll},legend_cell{l,ll})
    end
    end
    
end

%% Histograms of frequency band powers
% histogram plot settings: bw bin width, bd distance between groups of
% bins, bd distance in between bins.
bw=1;
bgd=2;
bbd=.2;

% figure for each Area
% subplots with each genotype-band for both absolute-normalized
uniq.Area_full=unique(bandtable.Area);
nArea_full=length(uniq.Area);
for ii=1:nArea_full
    fig.Area(ii)=figure('Name',['Band powers for ',uniq.Area_full{ii}]);
end

Genotype_bandtable=unique(bandtable.Genotype);
nGenotype_bandtable=length(Genotype_bandtable);

uniq.State=unique(bandtable.State);
nState=length(uniq.State);
power_cellmatrix=cell(nState,nArea_full,nGenotype_bandtable,nBand,2);
for i=1:nState
logical.i=strcmp(bandtable.State,uniq.State(i));
uniq.Area=unique(bandtable.Area(logical.i));
nArea=length(uniq.Area);
for ii=1:nArea
set(0,'CurrentFigure',fig.Area(ii))
logical.ii=logical.i&strcmp(bandtable.Area,uniq.Area(ii));
uniq.Genotype=unique(bandtable.Genotype(logical.ii));
% uniq.Genotype=["wt","het"]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nGenotype=length(uniq.Genotype);

%%% define bin location
temp.step2=bbd+bw;
temp.step1=bw+(nGenotype-1)*(temp.step2)+bgd;
temp.points=(bgd:temp.step1:bgd+(nBand-1)*temp.step1)';
bin_matrix=zeros(nBand,nGenotype);
for j=1:nBand
    bin_matrix(j,:)=temp.points(j)+bw/2:temp.step2:temp.points(j)+bw/2+(nGenotype-1)*temp.step2;
end

for g=1:nGenotype
logical.ii_g=logical.ii&strcmp(bandtable.Genotype,uniq.Genotype(g));
for jj=1:2
    
%%% start plot
subplot(nState,2,(i-1)*nState+jj), hold on
for j=1:nBand
    if jj==1
        temp.power=bandtable.Power{logical.ii_g}(j,:);
    else %normalize
        temp.power=bandtable.Power{logical.ii_g}(j,:)./bandtable.Power{logical.ii_g}(end,:);
    end
    ind.ii=find(strcmp(uniq.Area_full,uniq.Area(ii)));
    ind.g=find(strcmp(uniq.Genotype(g),Genotype_bandtable));
    power_cellmatrix{i,ind.ii,ind.g,j,jj}=temp.power;
    temp.power_mean=mean(temp.power);
    temp.SEM=std(temp.power)/sqrt(length(temp.power));
    bar(bin_matrix(j,g),temp.power_mean,...
    'FaceColor',c(g,:),'BarWidth',bw);
    errorbar(bin_matrix(j,g),temp.power_mean,temp.SEM,...
    'Color',c(g,:),'LineStyle','none','LineWidth',lw)
end

%xticks
xticks(temp.points+(bw+(nGenotype-1)*temp.step2)/2)
xticklabels(bandstr)

%legend
for gg=1:nGenotype
    h_bin(gg)=bar(NaN,'FaceColor',c(gg,:),'BarWidth',bw);
end
legend(h_bin,uniq.Genotype)

%title and ylabel
if jj==1
title([uniq.Area{ii},' ',uniq.State{i},' -absolute'])
ylabel('Power (\muV^2)')
else
title([uniq.Area{ii},' ',uniq.State{i},' -normalized'])
ylabel('Power (%)')
end
%%%% end plot

end
end
end
end

%Save statistics to table
count=0;
for i=1:nState 
for ii=1:nArea_full
for j=1:nBand
for jj=1:2
    
    if jj==1
        normstr="Absolute";
    else
        normstr="Normalized";
    end
    
H=cell(nGenotype_bandtable); 
P=H;
CI=H;
STATS=H;
for g=1:nGenotype_bandtable
for gg=1:nGenotype_bandtable
    
    if ~isempty(power_cellmatrix{i,ii,g,j,jj})&&~isempty(power_cellmatrix{i,ii,gg,j,jj})
    [H{g,gg},P{g,gg},CI{g,gg},STATS{g,gg}]...
    =ttest2(power_cellmatrix{i,ii,g,j,jj},power_cellmatrix{i,ii,gg,j,jj},...
    'vartype','unequal');
    end
    
end
end

    count=count+1;
    
    bandtable2.state(count,1)=uniq.State(i);
    bandtable2.area(count,1)=uniq.Area_full(ii);
    bandtable2.normalization(count,1)=normstr;
    bandtable2.band(count,1)=bandstr(j);
    
    temp.charH='H(:,1)';
    temp.charP='P(:,1)';
    temp.charCI='CI(:,1)';
    temp.charSTATS='STATS(:,1)';
    for ggg=2:nGenotype_bandtable
        temp.charH=[temp.charH,',H(:,',num2str(ggg),')'];
        temp.charP=[temp.charP,',P(:,',num2str(ggg),')'];
        temp.charCI=[temp.charCI,',CI(:,',num2str(ggg),')'];
        temp.charSTATS=[temp.charSTATS,',STATS(:,',num2str(ggg),')'];
    end
    bandtable2.H{count,1}=eval(['table(',temp.charH,')']);
    bandtable2.H{count,1}.Properties.RowNames={Genotype_bandtable{:}};
    bandtable2.H{count,1}.Properties.VariableNames={Genotype_bandtable{:}};
    bandtable2.P{count,1}=eval(['table(',temp.charP,')']);
    bandtable2.P{count,1}.Properties.RowNames={Genotype_bandtable{:}};
    bandtable2.P{count,1}.Properties.VariableNames={Genotype_bandtable{:}};
    bandtable2.CI{count,1}=eval(['table(',temp.charCI,')']);
    bandtable2.CI{count,1}.Properties.RowNames={Genotype_bandtable{:}};
    bandtable2.CI{count,1}.Properties.VariableNames={Genotype_bandtable{:}};
    bandtable2.STATS{count,1}=eval(['table(',temp.charSTATS,')']);
    bandtable2.STATS{count,1}.Properties.RowNames={Genotype_bandtable{:}};
    bandtable2.STATS{count,1}.Properties.VariableNames={Genotype_bandtable{:}};
    bandtable2.P_matrix{count,1}=cell2mat(P);
    
end
end
end
end
TableNames={'State','Area','Normalization','Band',...
            'Reject_Null_Hypothesis','p_value','CI','STATS','p_value_Matrix'};
BandTable=table(bandtable2.state,bandtable2.area,bandtable2.normalization,bandtable2.band,...
                bandtable2.H,bandtable2.P,bandtable2.CI,bandtable2.STATS,bandtable2.P_matrix,...
                'VariableNames',TableNames);
save('BandTable.mat', 'BandTable');