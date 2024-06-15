% Plot mean PSD for each mouse and genotype
clearvars, close all
%% Import
EF = uigetdir('','Select Experiment Folder (EF)');
cd(EF)
E=dir('*_*_*.mat');
nE=length(E);

load('DataTable')

% DataTable_cell1=cell(nE,1);
% DataTable_cell2=cell(nE,1);
% DataTable_cell3=cell(nE,1);
% E_FileName=cell(nE,1);
% for i=1:nE
%     E_FileName{i}=E(i).name(1:end-4);
%     load(E_FileName{i},'DataTable','DataTable_Current','DataTable_Current')
%     nM=size(DataTable,1);
%     checkifempty=true(nM,1);
%     for ii=1:nM
%         checkifempty(ii)=~isempty(DataTable.PSD{ii});
%     end
%     
%     DataTable_cell1{i}=DataTable(checkifempty,:);
%     DataTable_cell2{i}=DataTable_Current;
%     DataTable_cell3{i}=DataTable_Current;
% end
% DataTable1=cat(1,DataTable_cell1{:});
% DataTable_Current=cat(1,DataTable_cell2{:});
% DataTable_Current=cat(1,DataTable_cell3{:});

% Load and alter settings
load('Settings')
Settings.BandRanges(1,:)=[2,4];

% Assign current DataTable for plotting
DataTable_Current=DataTable;

% Reject CAMK mice-channels
LeaveOut=(...
    ...
    (strcmp(DataTable_Current.Name,'18210-03'))|... % reject all data
    (strcmp(DataTable_Current.Name,'18151-01'))|... % reject all data
    (strcmp(DataTable_Current.Name,'20363-03'))|... % reject all data
    (strcmp(DataTable_Current.Name,'Unknown1'))|... % reject all data
    ...
    (strcmp(DataTable_Current.Name,'16725-00')&strcmp(DataTable_Current.Day,'2011-11-09-T1-MaybeNoData'))|... % reject day(s)
    (strcmp(DataTable_Current.Name,'18113-02')&(strcmp(DataTable_Current.Day,'2011-09-14-T1-NotCertain')|strcmp(DataTable_Current.Day,'2011-09-15-T1')|strcmp(DataTable_Current.Day,'2011-09-23-T2')|strcmp(DataTable_Current.Day,'2011-09-27-T1')))|... % reject day(s)
    ...
    (strcmp(DataTable_Current.Name,'18113-00')&(strcmp(DataTable_Current.Day,'2011-10-19-T1')|strcmp(DataTable_Current.Day,'2011-09-13-T1-NotCertain')|strcmp(DataTable_Current.Day,'2011-10-03-T1')|strcmp(DataTable_Current.Day,'2011-10-07-T1'))&(strcmp(DataTable_Current.Channel,'left-M')|strcmp(DataTable_Current.Channel,'left-S')))|... % reject channel(s)
    (strcmp(DataTable_Current.Name,'18113-00')&~(strcmp(DataTable_Current.Day,'2011-10-03-T1')|strcmp(DataTable_Current.Day,'2011-10-07-T1')|strcmp(DataTable_Current.Day,'2011-09-15-T1'))&(strcmp(DataTable_Current.Channel,'right-M')|strcmp(DataTable_Current.Channel,'right-S')))|... % reject channel(s)
    (strcmp(DataTable_Current.Name,'20363-00')&(strcmp(DataTable_Current.Day,'2012-01-25-T1'))&(strcmp(DataTable_Current.Channel,'left-M')|strcmp(DataTable_Current.Channel,'left-S'))))... % reject channel(s)
    ;
DataTable_Current(LeaveOut,:)=[];

% note
% 15056-05_unknown_D1_Analyzed_noREM

p_threshold=.05;

%% Constants defined here (can be changed to input prompt)
maxFreq=50;
SEM_transparency=0.2;
normFreq=50;

% define fff by assuming fs to be constant
fs=500.2748;
fff=0:fs/round(fs*Settings.window_FFT):fs/2;

%% Plot settings
% plot property 'Color'
c=colormap(lines); % colors to most contrasting color map
% switch around colors (for wt and het)
% c(1:2,:)=c([2,1],:);%%%%%%%%%%%%%%%
c(1:3,:)=c([2,3,1],:);%%%%%%%%%%%%%%%
% c(2:3,:)=c([3,2],:);
close(gcf)

% plot property 'LineStyle'
ls{1}='-';
ls{2}='--';
ls{3}=':';
ls{4}='-.';

% plot property 'LineWidth'
lw=1.5;

%% Name each Day plots for each state and each channel
plot_logarithmic=2;
plot_PSDeachday=false;
maxFreq=20;
if plot_PSDeachday

uniq.Day_Full=unique(DataTable.Day);
uniq.State=unique(DataTable_Current.State);
for i=1:length(uniq.State)
temp.Logical_i=strcmp(DataTable_Current.State,uniq.State{i});
uniq.Channel=unique(DataTable_Current.Channel(temp.Logical_i));
    for ii=1:length(uniq.Channel)
    temp.Logical_ii=strcmp(DataTable_Current.Channel,uniq.Channel{ii})&temp.Logical_i;
    uniq.Name=unique(DataTable_Current.Name(temp.Logical_ii));
    
        for iii=1:length(uniq.Name)
        temp.Logical_iii=strcmp(DataTable_Current.Name,uniq.Name{iii})&temp.Logical_ii;
        uniq.Day=DataTable_Current.Day(temp.Logical_iii);
        figure('Name',[uniq.State{i},' ',uniq.Channel{ii},...
               ' Each Day for Name'])
        for j=1:plot_logarithmic
        subplot(1,2,j), hold on, grid on
        if j==1
            title('Absolute')
        else
            title('Relative')
        end
    
        temp.Power=DataTable_Current.PSD(temp.Logical_iii);
        
        % Plot if
        if sum(temp.Logical_iii)~=0
    
        hhh=zeros(sum(temp.Logical_iii),1);
        for iiii=1:sum(temp.Logical_iii)
            
        % for normalized plots
        if (j==2||j==4)
            temp.Power{iiii}=(temp.Power{iiii}/bandpower2(temp.Power{iiii},fff,[Settings.HP_ele,normFreq]))*100;
        end
        
        c_index=find(strcmp(uniq.Day_Full,uniq.Day{iiii}));
        while c_index>size(c,1)
            c_index=c_index-size(c,1);
        end
        ls_index=find(strcmp(uniq.Day_Full,uniq.Day{iiii}));
        while ls_index>length(ls)
            ls_index=ls_index-length(ls);
        end
        
        plot(fff,temp.Power{iiii},'Color',c(c_index,:),'LineStyle',ls{ls_index},'LineWidth',lw)
        hhh(iiii)=plot(NaN,NaN,'Color',c(c_index,:),'LineStyle',ls{ls_index},'LineWidth',lw);
            
        end
    
        % Figure settings
        legend(hhh, {uniq.Day{:}},'Interpreter','none');
        xlim([Settings.HP_ele,maxFreq])
%         xticks(linspace(0,maxFreq,6))
        title(uniq.Name{iii})
        xlabel('Frequency (Hz)')
        if j==3||j==4 % for logarithmic plots
            set(gca, 'YScale', 'log')
        end
        if j==1||j==3 % for non-normalized plots
            ylabel('LFPs PSD (\muV^2)')
        else % for normalized plots
            ylabel('LFPs PSD (%)')
        end
        
    
        end
        end % end of 1 plot
        end % end of 4 plots
    
    end
end

end

%% Plots (beside each day)

nState_full=length(unique(DataTable_Current.State));

uniq.Name_full=unique(DataTable_Current.Name);
uniq.State_full=unique(DataTable_Current.State);

count=0;
uniq.Channel=unique(DataTable_Current.Channel);
nChannel=length(uniq.Channel);
for i=1:nChannel
    logical.channel=strcmp(DataTable_Current.Channel,uniq.Channel(i));
    
    figure('Name',[uniq.Channel{i},' Mean over Days'])
    legend_cell=cell(nState_full,2);
    handle_cell=cell(nState_full,2);
    
    uniq.Name=unique(DataTable_Current.Name(logical.channel));
    
    for ii=1:length(uniq.Name)
    logical.channel_name=logical.channel&strcmp(DataTable_Current.Name,uniq.Name{ii});
    uniq.State=unique(DataTable_Current.State(logical.channel_name));
    
    for iii=1:length(uniq.State)
    logical.channel_name_state=logical.channel_name&strcmp(DataTable_Current.State,uniq.State{iii});
        
    if sum(logical.channel_name_state)
    
    count=count+1;    
        
    % do mean over days and save in cell for temp table
    PSD_Matrix=[DataTable_Current.PSD{logical.channel_name_state}];
    [temptable.PSD{count,1},SEM_meanday]=meanSEM(PSD_Matrix);
    
    % additional strings for temp table
    temptable.Channel(count,1)=uniq.Channel(i);
    temptable.Genotype(count,1)=DataTable_Current.Genotype(find(logical.channel_name_state,1,'first'));
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
    subplot(nState_full,2,(state_index-1)*2+1), grid on % index for plot (2x2 : state x norm)
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
    title([uniq.Channel{i},' ',uniq.State_full{state_index},' -absolute'])
    xlim([Settings.HP_ele,maxFreq])
    legend_cell{iii,1}=[legend_cell{iii,1},DataTable_Current.Name(find(logical.channel_name_state,1,'first'))];
    handle_cell{iii,1}=[handle_cell{iii,1},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    
    % normalization (%)
    PSD_Matrix=(PSD_Matrix./bandpower2(PSD_Matrix,fff,[Settings.HP_ele,normFreq]))*100; 
    [PSD_meanday_normalized,SEM_meanday_normalized]=meanSEM(PSD_Matrix);
    
    % normalized plot
    subplot(nState_full,2,(state_index-1)*2+2), grid on % index for plot
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
    title([uniq.Channel{i},' ',uniq.State_full{state_index},' -normalized'])
    xlim([Settings.HP_ele,maxFreq])
    legend_cell{iii,2}=[legend_cell{iii,2},DataTable_Current.Name(find(logical.channel_name_state,1,'first'))];
    handle_cell{iii,2}=[handle_cell{iii,2},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (%)')
    
    end % end of if-statement
    
    end
    end

for l=1:nState_full
    for ll=1:2
        subplot(nState_full,2,ll+(l-1)*2), grid on
        legend(handle_cell{l,ll},legend_cell{l,ll})
    end
end

end 

% clearvars handle_cell legend_cell
%% Genotype over Names
% in and after this section fs is assumed to be constant

% pre
count=0;

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
    subplot(nState_full,2,(state_index-1)*2+1), grid on
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
    
% temporary 'table'
count=count+1;
bandtable.Channel(count,1)=uniq.Channel(i);
bandtable.State(count,1)=uniq.State(ii);
bandtable.Genotype(count,1)=uniq.Genotype(iii);
for bbb=1:Settings.nBand
bandtable.PSD{count,1}(bbb,:)=bandpower2(PSD_Matrix,fff,Settings.BandRanges(bbb,:));
end

    % normalization (%)
    PSD_Matrix=(PSD_Matrix./bandpower2(PSD_Matrix,fff,[Settings.HP_ele,normFreq]))*100; 
    [temp.power_norm,temp.SEM_norm]=meanSEM(PSD_Matrix);
    
    % normalized plot
    subplot(nState_full,2,state_index*2), grid on
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
    end

for l=1:nState_full
    for ll=1:2
        subplot(nState_full,2,ll+(l-1)*2), grid on
        legend(handle_cell{l,ll},legend_cell{l,ll})
    end
end
    
end

%% Histograms of frequency band powers (do this for channels instead of channels)
% histogram plot settings: bw bin width, bd distance between groups of
% bins, bd distance in between bins.
bw=1;
bgd=2;
bbd=.2;

% figure for each Channel
% subplots with each genotype-band for both absolute-normalized
uniq.Channel_full=unique(bandtable.Channel);
nChannel_full=length(uniq.Channel);
for ii=1:nChannel_full
    fig.Channel(ii)=figure('Name',['Band powers for ',uniq.Channel_full{ii}]);
end

Genotype_bandtable=unique(bandtable.Genotype);
nGenotype_bandtable=length(Genotype_bandtable);

uniq.State=unique(bandtable.State);
nState=length(uniq.State);
power_cellmatrix=cell(nState,nChannel_full,nGenotype_bandtable,Settings.nBand,2);
for i=1:nState
logical.i=strcmp(bandtable.State,uniq.State(i));
uniq.Channel=unique(bandtable.Channel(logical.i));
nChannel=length(uniq.Channel);
for ii=1:nChannel
set(0,'CurrentFigure',fig.Channel(ii))
logical.ii=logical.i&strcmp(bandtable.Channel,uniq.Channel(ii));
uniq.Genotype=unique(bandtable.Genotype(logical.ii));
% uniq.Genotype=["wt","het"]; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nGenotype=length(uniq.Genotype);

%%% define bin location
temp.step2=bbd+bw;
temp.step1=bw+(nGenotype-1)*(temp.step2)+bgd;
temp.points=(bgd:temp.step1:bgd+(Settings.nBand-1)*temp.step1)';
bin_matrix=zeros(Settings.nBand,nGenotype);
for j=1:Settings.nBand
    bin_matrix(j,:)=temp.points(j)+bw/2:temp.step2:temp.points(j)+bw/2+(nGenotype-1)*temp.step2;
end

for g=1:nGenotype
logical.ii_g=logical.ii&strcmp(bandtable.Genotype,uniq.Genotype(g));
for jj=1:2
    
%%% start plot
subplot(nState,2,(i-1)*2+jj), hold on, grid on
for j=1:Settings.nBand
    if jj==1
        temp.power=bandtable.PSD{logical.ii_g}(j,:);
    else %normalize
        temp.power=bandtable.PSD{logical.ii_g}(j,:)./bandtable.PSD{logical.ii_g}(end,:);
    end
    ind.ii=find(strcmp(uniq.Channel_full,uniq.Channel(ii)));
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
xticklabels(Settings.Bands)

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
for j=1:Settings.nBand
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
    bandtable2.channel(count,1)=uniq.Channel_full(ii);
    bandtable2.normalization(count,1)=normstr;
    bandtable2.band(count,1)=Settings.Bands(j);
    
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
    bandtable2.Significance(count,1)=round(sum(sum(cell2mat(P)<p_threshold))/2);
    
end
end
end
end
TableNames={'State','Channel','Normalization','Band',...
            'Reject_Null_Hypothesis','p_value','CI','STATS','p_value_Matrix','Significance'};
BandTable=table(bandtable2.state,bandtable2.channel,bandtable2.normalization,bandtable2.band,...
                bandtable2.H,bandtable2.P,bandtable2.CI,bandtable2.STATS,bandtable2.P_matrix,bandtable2.Significance,...
                'VariableNames',TableNames);
save('BandTable.mat', 'BandTable');