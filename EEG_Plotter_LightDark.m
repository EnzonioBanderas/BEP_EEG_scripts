% Plot mean PSD for each mouse and genotype
clearvars, close all
%% Import
EF = uigetdir('','Select Experiment Folder (EF)');
cd(EF)
E=dir('*_Analyzed.mat');
nE=length(E);
load('DataTable')

load('Settings')
DataTable5=1;
% Leave out high cross correlation between channels
LeaveOut=strcmp(DataTable2.Name,'14248-05')&(strcmp(DataTable2.Channel,'right-M')|strcmp(DataTable2.Channel,'right-S'))|...
        strcmp(DataTable2.Name,'14346-01')&(strcmp(DataTable2.Channel,'left-S') |strcmp(DataTable2.Channel,'left-M')) |...
        strcmp(DataTable2.Name,'14346-03')&(strcmp(DataTable2.Channel,'right-M')|strcmp(DataTable2.Channel,'right-S')|strcmp(DataTable2.Channel,'right-S') |strcmp(DataTable2.Channel,'right-M'))|...
        strcmp(DataTable2.Name,'14351-02')&(strcmp(DataTable2.Channel,'right-S')|strcmp(DataTable2.Channel,'right-M'))|...
        strcmp(DataTable2.Name,'14463-02')&(strcmp(DataTable2.Channel,'right-M')|strcmp(DataTable2.Channel,'right-S'));
DataTable2(LeaveOut,:)=[];

%% Constants defined here (can be changed to input prompt)
% Define between which hours it is light
LightDark=[7,19];

maxFreq=20;
SEM_transparency=0.2;
normFreq=50;

% define fff by assuming fs to be constant
load(E(1).name,'header')
fs=header.frequency(1);
fff=0:fs/round(fs*Settings.window_FFT):fs/2;

%% Plot settings
% plot property 'Color'
c=colormap(lines); % colors to most contrasting color map
% switch around colors (for wt and het)
c(1:3,:)=c([2,3,1],:);%%%%%%%%%%%%%%%
close(gcf)

% plot property 'LineStyle'
ls{1}='-';
ls{2}='--';
ls{3}=':';
ls{4}='-.';

% plot property 'LineWidth'
lw=1.5;

Bands=[Settings.Bands;{'Full_Kreuzer'};{'Full'}];
nBand=length(Bands);
States={'Awake_Light';'NREM_Light';'REM_Light';...
        'Awake_Dark';'NREM_Dark';'REM_Dark'};
nState=length(States);

%% Plot and create temptable

nState_full=length(unique(DataTable2.State));

uniq.Name_full=unique(DataTable2.Name);
uniq.State_full=unique(DataTable2.State);

count=0;
uniq.Channel=unique(DataTable2.Channel);
nChannel=length(uniq.Channel);
for i=1:nChannel
    logical.channel=strcmp(DataTable2.Channel,uniq.Channel(i));
    
    figure('Name',[uniq.Channel{i},' Mean over Days'])
    legend_cell=cell(nState_full,2);
    handle_cell=cell(nState_full,2);
    
    uniq.Name=unique(DataTable2.Name(logical.channel));
    
    for ii=1:length(uniq.Name)
    logical.channel_name=logical.channel&strcmp(DataTable2.Name,uniq.Name{ii});
    uniq.State=unique(DataTable2.State(logical.channel_name));
    
    for iii=1:length(uniq.State)
    logical.channel_name_state=logical.channel_name&strcmp(DataTable2.State,uniq.State{iii});
        
    if sum(logical.channel_name_state)
    
    count=count+1;    
        
    % do mean over days and save in cell for temp table
    PSD_Matrix=[DataTable2.PSD{logical.channel_name_state}];
    [temptable.PSD{count,1},SEM_meanday]=meanSEM(PSD_Matrix);
    
    % additional strings for temp table
    temptable.Channel(count,1)=uniq.Channel(i);
    temptable.Genotype(count,1)=DataTable2.Genotype(find(logical.channel_name_state,1,'first'));
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
    subplot(nState,2,(state_index-1)*2+1) % index for plot (2x2 : state x norm)
    h=plot(fff,temptable.PSD{count},...
        'Color',c(name_index_c,:),'LineStyle',ls{name_index_ls},'LineWidth',lw); hold on
    patch([fff,...
           fff(end:-1:1),...
           fff(1)],...
          [temptable.PSD{count}-SEM_meanday;...
           temptable.PSD{count}(end:-1:1)+SEM_meanday(end:-1:1);...
           temptable.PSD{count}(1)],...
           c(name_index_c,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
    title([uniq.Channel{i},' ',uniq.State_full{state_index},' -absolute'],'Interpreter','none')
    xlim([Settings.HP_ele,maxFreq])
    legend_cell{iii,1}=[legend_cell{iii,1},DataTable2.Name(find(logical.channel_name_state,1,'first'))];
    handle_cell{iii,1}=[handle_cell{iii,1},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    
    % normalization (%)
    PSD_Matrix=(PSD_Matrix./bandpower2(PSD_Matrix,fff,[Settings.HP_ele,normFreq]))*100; 
    [PSD_meanday_normalized,SEM_meanday_normalized]=meanSEM(PSD_Matrix);
    
    % normalized plot
    subplot(nState,2,state_index*2) % index for plot
    h=plot(fff,PSD_meanday_normalized,...
        'Color',c(name_index_c,:),'LineStyle',ls{name_index_ls},'LineWidth',lw); hold on
    patch([fff,...
           fff(end:-1:1),...
           fff(1)],...
          [PSD_meanday_normalized-SEM_meanday_normalized;...
           PSD_meanday_normalized(end:-1:1)+SEM_meanday_normalized(end:-1:1);...
           PSD_meanday_normalized(1)],...
           c(name_index_c,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency)
    title([uniq.Channel{i},' ',uniq.State_full{state_index},' -normalized'],'Interpreter','none')
    xlim([Settings.HP_ele,maxFreq])
    legend_cell{iii,2}=[legend_cell{iii,2},DataTable2.Name(find(logical.channel_name_state,1,'first'))];
    handle_cell{iii,2}=[handle_cell{iii,2},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (%)')
    
    end % end of if-statement
    
    end
    end

% for l=1:nState_full
%     for ll=1:2
%         subplot(nState_full,2,ll+(l-1)*2)
%         legend(handle_cell{l,ll},legend_cell{l,ll})
%     end
% end
        subplot(nState_full,2,1)
        LEG=legend(handle_cell{1,1},legend_cell{1,1});
        LEG.Position=[0.4807 0.4929 0.0608 0.0578];

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
    PSD_Matrix=cat(2,temptable.PSD{logical.channel_state_genotype});
    [temp.power,temp.SEM]=meanSEM(PSD_Matrix);
    
    if ~isempty(temp.power)
        
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
    title([uniq.Channel{i},' ',uniq.State_full{state_index},' -absolute'],'Interpreter','none')
    xlim([Settings.HP_ele,maxFreq])
    legend_cell{ii,1}=[legend_cell{ii,1},temptable.Genotype(find(logical.channel_state_genotype,1,'first'))];
    handle_cell{ii,1}=[handle_cell{ii,1},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    
    % normalization (%)
    PSD_Matrix=(PSD_Matrix./bandpower2(PSD_Matrix,fff,[Settings.HP_ele,normFreq]))*100; 
    [temp.power_norm,temp.SEM_norm]=meanSEM(PSD_Matrix);
    
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
    title([uniq.Channel{i},' ',uniq.State_full{state_index},' -normalized'],'Interpreter','none')
    xlim([Settings.HP_ele,maxFreq])
    legend_cell{ii,2}=[legend_cell{ii,2},temptable.Genotype(find(logical.channel_state_genotype,1,'first'))];
    handle_cell{ii,2}=[handle_cell{ii,2},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (%)')
    
    end
    
    end
    end
    
    subplot(nState_full,2,1)
    LEG=legend(handle_cell{1},legend_cell{1});
    LEG.Position=[0.4807 0.4929 0.0608 0.0578];
    
end
      

%% Histograms of frequency band powers (do this for channels instead of areas)
% histogram plot settings: bw bin width, bd distance between groups of
% bins, bd distance in between bins.
bw=1;
bgd=2;
bbd=.2;

% figure for each Channel
% subplots with each genotype-band for both absolute-normalized
uniq.Channel_full=unique(temptable.Channel);
nChannel_full=length(uniq.Channel);
for ii=1:nChannel_full
    fig.Channel(ii)=figure('Name',['Band powers for ',uniq.Channel_full{ii}]);
end

uniq.Genotype_full=unique(temptable.Genotype);
nGenotype_full=length(uniq.Genotype_full);

uniq.State=unique(temptable.State);
nState=length(uniq.State);
power_cellmatrix=cell(nState,nChannel_full,nGenotype,nBand,2);
for i=1:nState
logical.i=strcmp(temptable.State,uniq.State(i));
uniq.Channel=unique(temptable.Channel(logical.i));
nChannel=length(uniq.Channel);
for ii=1:nChannel
set(0,'CurrentFigure',fig.Channel(ii))
logical.ii=logical.i&strcmp(temptable.Channel,uniq.Channel(ii));
uniq.Genotype=unique(temptable.Genotype(logical.ii));
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
logical.ii_g=logical.ii&strcmp(temptable.Genotype,uniq.Genotype(g));
for jj=1:2
    
%%% start plot
subplot(nState,2,(i-1)*2+jj), hold on
for j=1:nBand
    if jj==1
        temp.power=temptable.PSD{logical.ii_g}(j,:);
    else %normalize
        temp.power=temptable.PSD{logical.ii_g}(j,:)./temptable.PSD{logical.ii_g}(end,:);
    end
    ind.ii=find(strcmp(uniq.Channel_full,uniq.Channel(ii)));
    ind.g=find(strcmp(uniq.Genotype(g),uniq.Genotype));
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
xticklabels(Bands)

%legend
for gg=1:nGenotype
    h_bin(gg)=bar(NaN,'FaceColor',c(gg,:),'BarWidth',bw);
end
legend(h_bin,uniq.Genotype,'Interpreter','none')

%title and ylabel
if jj==1
title([uniq.Channel{ii},' ',uniq.State{i},' -absolute'])
ylabel('PSD (\muV^2)')
else
title([uniq.Channel{ii},' ',uniq.State{i},' -normalized'])
ylabel('PSD (%)')
end
%%%% end plot

end
end
end
end

%Save statistics to table
count=0;
for i=1:nState 
for ii=1:nChannel_full
for j=1:nBand
for jj=1:2
    
    if jj==1
        normstr="Absolute";
    else
        normstr="Normalized";
    end
    
H=cell(nuniq.Genotype); 
P=H;
CI=H;
STATS=H;
for g=1:nuniq.Genotype
for gg=1:nuniq.Genotype
    
    if ~isempty(power_cellmatrix{i,ii,g,j,jj})&&~isempty(power_cellmatrix{i,ii,gg,j,jj})
    [H{g,gg},P{g,gg},CI{g,gg},STATS{g,gg}]...
    =ttest2(power_cellmatrix{i,ii,g,j,jj},power_cellmatrix{i,ii,gg,j,jj},...
    'vartype','unequal');
    end
    
end
end

    count=count+1;
    
    temptable2.state(count,1)=uniq.State(i);
    temptable2.area(count,1)=uniq.Channel_full(ii);
    temptable2.normalization(count,1)=normstr;
    temptable2.band(count,1)=Bands(j);
    
    temp.charH='H(:,1)';
    temp.charP='P(:,1)';
    temp.charCI='CI(:,1)';
    temp.charSTATS='STATS(:,1)';
    for ggg=2:nuniq.Genotype
        temp.charH=[temp.charH,',H(:,',num2str(ggg),')'];
        temp.charP=[temp.charP,',P(:,',num2str(ggg),')'];
        temp.charCI=[temp.charCI,',CI(:,',num2str(ggg),')'];
        temp.charSTATS=[temp.charSTATS,',STATS(:,',num2str(ggg),')'];
    end
    temptable2.H{count,1}=eval(['table(',temp.charH,')']);
    temptable2.H{count,1}.Properties.RowNames={uniq.Genotype{:}};
    temptable2.H{count,1}.Properties.VariableNames={uniq.Genotype{:}};
    temptable2.P{count,1}=eval(['table(',temp.charP,')']);
    temptable2.P{count,1}.Properties.RowNames={uniq.Genotype{:}};
    temptable2.P{count,1}.Properties.VariableNames={uniq.Genotype{:}};
    temptable2.CI{count,1}=eval(['table(',temp.charCI,')']);
    temptable2.CI{count,1}.Properties.RowNames={uniq.Genotype{:}};
    temptable2.CI{count,1}.Properties.VariableNames={uniq.Genotype{:}};
    temptable2.STATS{count,1}=eval(['table(',temp.charSTATS,')']);
    temptable2.STATS{count,1}.Properties.RowNames={uniq.Genotype{:}};
    temptable2.STATS{count,1}.Properties.VariableNames={uniq.Genotype{:}};
    temptable2.P_matrix{count,1}=cell2mat(P);
    
end
end
end
end
TableNames={'State','Channel','Normalization','Band',...
            'Reject_Null_Hypothesis','p_value','CI','STATS','p_value_Matrix'};
BandTable=table(temptable2.state,temptable2.area,temptable2.normalization,temptable2.band,...
                temptable2.H,temptable2.P,temptable2.CI,temptable2.STATS,temptable2.P_matrix,...
                'VariableNames',TableNames);
save('BandTable.mat', 'BandTable');



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

%%
% for i=1:nGenotype
% ttest2('vartype','unequal')