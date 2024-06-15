% Plot mean PSD for each mouse and genotype
clearvars, close all
%% Import
EF = uigetdir('','Select Experiment Folder (EF)');
cd(EF)
E=dir('*_Analyzed.mat');
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
Settings.BandRanges=Settings.BandRanges(1:end-2,:);
% Settings.Bands=Settings.Bands(1:end-2,:);
Settings.nBand=length(Settings.Bands);
% Settings.Bands=[Settings.Bands;{'Full_Kreuzer'};{'Full'}];
% Settings.nBand=length(Settings.Bands);

Settings.Bands=[Settings.Bands;{'Full'}];
Settings.BandRanges=[Settings.BandRanges;.5,100];
Settings.nBand=length(Settings.Bands);

% Assign current DataTable for plotting
DataTable_Current=DataTable5;

% Leave out high cross correlation between channels and no signal (less strict)
DataTable_Current.Genotype(strcmp(DataTable_Current.Genotype,'wtrescue'))={'wt'};%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LeaveOut=strcmp(DataTable_Current.Name,'14248-05')&(strcmp(DataTable_Current.Channel,'right-M')|strcmp(DataTable_Current.Channel,'right-S'))|...
         strcmp(DataTable_Current.Name,'14346-01')&(strcmp(DataTable_Current.Channel,'left-S') |strcmp(DataTable_Current.Channel,'left-M')) |...
         strcmp(DataTable_Current.Name,'14346-03')&(strcmp(DataTable_Current.Channel,'right-M')|strcmp(DataTable_Current.Channel,'right-S'))|...
         strcmp(DataTable_Current.Name,'14351-02')&(strcmp(DataTable_Current.Channel,'right-S')|strcmp(DataTable_Current.Channel,'right-M'))|...
         strcmp(DataTable_Current.Name,'14463-02')&(strcmp(DataTable_Current.Channel,'right-M')|strcmp(DataTable_Current.Channel,'right-S'))|...
         strcmp(DataTable_Current.Name,'14463-01')|... % both electrodes in both hemispheres have high CC within hemisphere
         strcmp(DataTable_Current.Name,'14463-03')&strcmp(DataTable_Current.Channel,'right-S')|... % no signal for this channel
         strcmp(DataTable_Current.Name,'14578-01')|... % Probably switched around electrodes
         strcmp(DataTable_Current.Name,'14783-02')&(strcmp(DataTable_Current.Channel,'right-S')|strcmp(DataTable_Current.Channel,'right-M'))|... % high CC
         strcmp(DataTable_Current.Name,'15096-05')&(strcmp(DataTable_Current.Channel,'left-S')|strcmp(DataTable_Current.Channel,'left-M'))|... % high CC
         strcmp(DataTable_Current.Name,'15096-05')&strcmp(DataTable_Current.Channel,'right-S')|...
         ... % G5
         strcmp(DataTable_Current.Name,'15332-05')|... % suggestion rejected everything (very high EEG amplitude everywhere)
         strcmp(DataTable_Current.Name,'15252-04')&(strcmp(DataTable_Current.Channel,'right-S')|strcmp(DataTable_Current.Channel,'right-M'))|... % white noise in right hemisphere
         strcmp(DataTable_Current.Name,'15250-05')|... % Decision 2 failed
         strcmp(DataTable_Current.Name,'15252-02')&strcmp(DataTable_Current.Channel,'right-S')|... % white noise in right-S
         ... % G6, where new standard where it is assumed that there are 'assymetric' copies is not enforced
         strcmp(DataTable_Current.Name,'16952-01')|...% special case where it seems that either M's are highly similar and one M is copied to right-S or right-S has been copied to both M's
         strcmp(DataTable_Current.Name,'16952-02')&(strcmp(DataTable_Current.Channel,'right-M')|strcmp(DataTable_Current.Channel,'right-S'))|... % high CC uncertain direction
         strcmp(DataTable_Current.Name,'16952-03')&(strcmp(DataTable_Current.Channel,'right-M')|strcmp(DataTable_Current.Channel,'right-S'))|... % high CC uncertain direction
         strcmp(DataTable_Current.Name,'16952-05')|... % special case where left-M and right-S seem to be one hemisphere pair and the other two seem to be another pair
         strcmp(DataTable_Current.Name,'16952-06')&~strcmp(DataTable_Current.Channel,'left-S')|... % no signal in rest of channels
         strcmp(DataTable_Current.Name,'16965-03')&(strcmp(DataTable_Current.Channel,'right-M')|strcmp(DataTable_Current.Channel,'right-S'))|... % no signal
         strcmp(DataTable_Current.Name,'17148-03')|...% special case where it seems that either M's are highly similar and one M is copied to left-S or left-S has been copied to both M's
     ...%strcmp(DataTable_Current.Name,'17148-06')&strcmp(DataTable_Current.Channel,'left-S')|... % no signal
         strcmp(DataTable_Current.Name,'17148-06')|... % no signal in left-S and overall high low frequency power
         strcmp(DataTable_Current.Name,'17928-03'); % not enough signal overall due to large amounts of noise
DataTable_Current(LeaveOut,:)=[];

% % reject all channels if there is any suspicion (very strict)
% % DataTable_Current.Genotype(strcmp(DataTable_Current.Genotype,'wtrescue'))={'wt'};%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LeaveOut=strcmp(DataTable_Current.Name,'14248-05')|...
%          strcmp(DataTable_Current.Name,'14346-01')|...
%          strcmp(DataTable_Current.Name,'14346-03')|...
%          strcmp(DataTable_Current.Name,'14351-02')|...
%          strcmp(DataTable_Current.Name,'14463-02')|...
%          strcmp(DataTable_Current.Name,'14463-01')|... % both electrodes in both hemispheres have high CC within hemisphere
%          strcmp(DataTable_Current.Name,'14463-03')|... % no signal for this channel
%          strcmp(DataTable_Current.Name,'14578-01')|... % Probably switched around electrodes
%          strcmp(DataTable_Current.Name,'14783-02')|... % high CC
%          strcmp(DataTable_Current.Name,'15056-05')|... % NoREM
%          strcmp(DataTable_Current.Name,'15056-01')|... % high CC M
%          strcmp(DataTable_Current.Name,'15096-04')|... % high CC M
%          strcmp(DataTable_Current.Name,'15096-05')|... % high CC and nosignal         
%          ... % G5
%          strcmp(DataTable_Current.Name,'15332-05')|... % suggestion rejected everything (very high EEG amplitude everywhere)
%          strcmp(DataTable_Current.Name,'15252-04')|... % white noise in right hemisphere
%          strcmp(DataTable_Current.Name,'15250-05')|... % Decision 2 failed
%          strcmp(DataTable_Current.Name,'15252-02')|... % white noise in right-S
%          strcmp(DataTable_Current.Name,'15252-05'); % suggestion probably removes most of trace
% DataTable_Current(LeaveOut,:)=[];

% note
% 15056-05_unknown_D1_Analyzed_noREM

p_threshold=.05;

%% Constants defined here (can be changed to input prompt)
maxFreq=50;
SEM_transparency=0.2;
normFreq=50;

% define fff by assuming fs to be constant
fs=204.8;
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

%% 

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

markerSize=5;
heightSig=1;

BarC= [     0    0.4470    0.7410;
       0.9290    0.6940    0.1250;
       0.8500    0.3250    0.0980;
       0.4940    0.1840    0.5560;
       0.4660    0.6740    0.1880;
       0.3010    0.7450    0.9330];


% figure for each Channel
% subplots with each genotype-band for both absolute-normalized
uniq.Channel_full=unique(bandtable.Channel);
nChannel_full=length(uniq.Channel);
for ii=1:nChannel_full
    fig.Channel(ii)=figure('Name',['Band powers for ',uniq.Channel_full{ii}]);
end

Genotype_bandtable=unique(bandtable.Genotype);
Genotype_bandtable={'wt';'hetrescue';'het'}; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change order
% uniq.Genotype={'wt';'wtrescue';'hetrescue';'het'}; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nGenotype_bandtable=length(Genotype_bandtable);

barHG=zeros(nGenotype_bandtable,1);

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
Genotype_bandtable={'wt';'hetrescue';'het'}; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% change order
% uniq.Genotype={'wt';'wtrescue';'hetrescue';'het'}; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        temp.power=(bandtable.PSD{logical.ii_g}(j,:)./bandtable.PSD{logical.ii_g}(end,:))*100;
    end
    ind.ii=find(strcmp(uniq.Channel_full,uniq.Channel(ii)));
    ind.g=find(strcmp(uniq.Genotype(g),Genotype_bandtable));
    power_cellmatrix{i,ind.ii,ind.g,j,jj}=temp.power;
    temp.power_mean=mean(temp.power);
    temp.SEM=std(temp.power)/sqrt(length(temp.power));
    BARS(ii,i,jj,j,g)=bar(bin_matrix(j,g),temp.power_mean,...
    'FaceColor',BarC(g,:),'BarWidth',bw);
    BARS_error(ii,i,jj,j,g)=errorbar(bin_matrix(j,g),temp.power_mean,temp.SEM,...
    'Color',BarC(g,:),'LineStyle','none','LineWidth',lw);
end

%xticks
xticks(temp.points+(bw+(nGenotype-1)*temp.step2)/2)
xticklabels(Settings.Bands)

%legend
for gg=1:nGenotype
    h_bin(gg)=bar(NaN,'FaceColor',BarC(gg,:),'BarWidth',bw);
end
legend(h_bin,uniq.Genotype,'Interpreter','none','AutoUpdate','off')

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
P=zeros(nGenotype_bandtable);
CI=H;
STATS=H;
check_matrix=triu(true(nGenotype),1);
for g=1:nGenotype_bandtable
for gg=1:nGenotype_bandtable
    
    if ~isempty(power_cellmatrix{i,ii,g,j,jj})&&~isempty(power_cellmatrix{i,ii,gg,j,jj})
    [H{g,gg},P(g,gg),CI{g,gg},STATS{g,gg}]...
    =ttest2(power_cellmatrix{i,ii,g,j,jj},power_cellmatrix{i,ii,gg,j,jj},...
    'vartype','unequal');
    end
    
end
end

% Pcheck=P<0.05&check_matrix;
% countSig=0;
% % fprintf('0')
% if any(Pcheck(:))
%     set(0,'CurrentFigure',fig.Channel(ii))
%     subplot(nState,2,(i-1)*2+jj)
% %     fprintf('F')
%     for g=1:nGenotype_bandtable
%         barHG(g)=BARS(ii,i,jj,j,g).YData;
%     end
%     barHG=max(barHG);
%     
%     for g=1:nGenotype_bandtable
%     for gg=1:nGenotype_bandtable
%     if Pcheck(g,gg)
% %         fprintf('S')
%         countSig=countSig+1;
%         barX=[BARS(ii,i,jj,j,g).XData,BARS(ii,i,jj,j,gg).XData];
% %         plot(1,1,'Marker','o')
%         barY=[BARS(ii,i,jj,j,g).YData,BARS(ii,i,jj,j,gg).YData];
%         heightSig2=countSig*heightSig*barHG/2;
%         errorY=[BARS_error(ii,i,jj,j,g).YPositiveDelta,BARS_error(ii,i,jj,j,gg).YPositiveDelta];
%         plot(barX,ones(1,2)*(max(barY)+heightSig2),'Color',[0.5,0.5,0.5],'LineWidth',2)%%
%         plot(barX(1)*ones(1,2),[barY(1)+errorY(1),max(barY)+heightSig2],'Color',[0.5,0.5,0.5],'LineWidth',2)
%         plot(barX(2)*ones(1,2),[barY(2)+errorY(2),max(barY)+heightSig2],'Color',[0.5,0.5,0.5],'LineWidth',2)
%         if P(g,gg)<0.05
%             sigStar(countSig)=plot(mean(barX),max(barY)+heightSig2,'Color',[0,0,0],'Marker','*','MarkerSize',markerSize);
%         elseif P(g,gg)<0.01
%             sigStar(countSig)=plot(mean(barX)-markerSize/2,max(barY)+heightSig2,'Color',[0,0,0],'Marker','s','MarkerSize',markerSize);
% %             plot(mean(barX)+markerSize/2,max(barY)+heightSig2,'Color',[0,0,0],'Marker','*','MarkerSize',markerSize)
%         elseif P(g,gg)<0.001
%             sigStar(countSig)=plot(mean(barX)-markerSize,max(barY)+heightSig2,'Color',[0,0,0],'Marker','o','MarkerSize',markerSize);
% %             plot(mean(barX),           max(barY)+heightSig2,'Color',[0,0,0],'Marker','*','MarkerSize',markerSize)
% %             plot(mean(barX)+markerSize,max(barY)+heightSig2,'Color',[0,0,0],'Marker','*','MarkerSize',markerSize)
%         end
%     end
%     end
%     end
%     for y=1:countSig
%         uistack(sigStar(y),'top')
%     end
%     clearvars sigStar
% end

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
    bandtable2.P_matrix{count,1}=P;
    bandtable2.Significance(count,1)=round(sum(sum(P<p_threshold))/2);
    
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

%% Analysis of state differences between genotypes %%
% 1) Percentages Awake NREM REM
% 2) Mean of distributions of lengths
% 3) Standard deviation of distributions

% Assign current DataTable for plotting
DataTable_Current=DataTable3;

% % Leave out high cross correlation between channels
% LeaveOut=strcmp(DataTable_Current.Name,'14248-05')|...
%          strcmp(DataTable_Current.Name,'14346-01')|...
%          strcmp(DataTable_Current.Name,'14346-03')|...
%          strcmp(DataTable_Current.Name,'14351-02')|...
%          strcmp(DataTable_Current.Name,'14463-02')|...
%          strcmp(DataTable_Current.Name,'14463-01')|... % both electrodes in both hemispheres have high CC within hemisphere
%          strcmp(DataTable_Current.Name,'14463-03')|... % no signal for this channel
%          strcmp(DataTable_Current.Name,'14578-01'); % Probably switched around electrodes
% DataTable_Current(LeaveOut,:)=[];

% % Leave out high cross correlation between channels
% LeaveOut=strcmp(DataTable_Current.Name,'12255-05');
% DataTable_Current(LeaveOut,:)=[];

% Leave out predictions without REM (where 2nd decision failed)
DataTable_Current.Genotype(strcmp(DataTable_Current.Genotype,'wtrescue'))={'wt'};%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LeaveOut=strcmp(DataTable_Current.Name,'15056-05')|... % 2nd decision fail
         strcmp(DataTable_Current.Name,'15250-05')|... % 2nd decision fail
         strcmp(DataTable_Current.Name,'16952-06'); % 2nd decision fail
DataTable_Current(LeaveOut,:)=[];

% Create DataMatrix from DataTable
count=0;
uniq.Name=unique(DataTable_Current.Name);
nName=length(uniq.Name);
clearvars Name_cell Genotype_cell State_cell LengthDistribution_cell
for i=1:nName
    logical.Name=strcmp(DataTable_Current.Name,uniq.Name{i});
    uniq.State=unique(DataTable_Current.State(logical.Name));
    nState=length(uniq.State);
    for ii=1:nState
        logical.Name_State=logical.Name&strcmp(DataTable_Current.State,uniq.State{ii});
        if any(logical.Name_State)
            count=count+1;
            Name_cell(count,1)=DataTable_Current.Name(find(logical.Name_State,1,'first'));
            Genotype_cell(count,1)=DataTable_Current.Genotype(find(logical.Name_State,1,'first'));
            State_cell(count,1)=DataTable_Current.State(find(logical.Name_State,1,'first'));
            LengthDistribution_cell{count,1}=cat(1,DataTable_Current.LengthDistribution{logical.Name_State});
        end
    end
end

DataMatrix=zeros(nName,18);
DataMatrix_Genotype=cell(nName,1);
for i=1:nName
    
    logical.Name=strcmp(Name_cell,uniq.Name{i});
    
    DataMatrix_Genotype(i)=Genotype_cell(find(logical.Name,1,'first'));
    
    DataMatrix(i,1)=sum(LengthDistribution_cell{logical.Name&strcmp(State_cell,'Awake_Light')});
    DataMatrix(i,2)=sum(LengthDistribution_cell{logical.Name&strcmp(State_cell,'NREM_Light')});
    DataMatrix(i,3)=sum(LengthDistribution_cell{logical.Name&strcmp(State_cell,'REM_Light')});
    DataMatrix(i,1:3)=(DataMatrix(i,1:3)/sum(DataMatrix(i,1:3)))*100;
    DataMatrix(i,4)=sum(LengthDistribution_cell{logical.Name&strcmp(State_cell,'Awake_Dark')});
    DataMatrix(i,5)=sum(LengthDistribution_cell{logical.Name&strcmp(State_cell,'NREM_Dark')});
    DataMatrix(i,6)=sum(LengthDistribution_cell{logical.Name&strcmp(State_cell,'REM_Dark')});
    DataMatrix(i,4:6)=(DataMatrix(i,4:6)/sum(DataMatrix(i,4:6)))*100;
    
    DataMatrix(i,7)=mean(LengthDistribution_cell{logical.Name&strcmp(State_cell,'Awake_Light')});
    DataMatrix(i,8)=mean(LengthDistribution_cell{logical.Name&strcmp(State_cell,'NREM_Light')});
    DataMatrix(i,9)=mean(LengthDistribution_cell{logical.Name&strcmp(State_cell,'REM_Light')});
    DataMatrix(i,10)=mean(LengthDistribution_cell{logical.Name&strcmp(State_cell,'Awake_Dark')});
    DataMatrix(i,11)=mean(LengthDistribution_cell{logical.Name&strcmp(State_cell,'NREM_Dark')});
    DataMatrix(i,12)=mean(LengthDistribution_cell{logical.Name&strcmp(State_cell,'REM_Dark')});
    
    DataMatrix(i,13)=std(LengthDistribution_cell{logical.Name&strcmp(State_cell,'Awake_Light')});
    DataMatrix(i,14)=std(LengthDistribution_cell{logical.Name&strcmp(State_cell,'NREM_Light')});
    DataMatrix(i,15)=std(LengthDistribution_cell{logical.Name&strcmp(State_cell,'REM_Light')});
    DataMatrix(i,16)=std(LengthDistribution_cell{logical.Name&strcmp(State_cell,'Awake_Dark')});
    DataMatrix(i,17)=std(LengthDistribution_cell{logical.Name&strcmp(State_cell,'NREM_Dark')});
    DataMatrix(i,18)=std(LengthDistribution_cell{logical.Name&strcmp(State_cell,'REM_Dark')});
    
end

% Plot all 18 histograms
H=cell(18,1);
P=cell(18,1);
P_Significance=zeros(18,1);
CI=cell(18,1);
STATS=cell(18,1);

% col=@(x)reshape(x,numel(x),1);
% boxplot2=@(C,varargin)boxplot(...
%          cell2mat(cellfun(col,col(C),'uni',0)),....
%          cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),...
%          varargin{:});

for i=1:18
    figure('Name',['DataMatrix Column ',num2str(i)])
    hold on
    grid on
%     h(1,1)=histogram(DataMatrix(strcmp(DataMatrix_Genotype,'wt'),i),'FaceColor',c(3,:),'NumBins',10);
%     h(2,1)=histogram(DataMatrix(strcmp(DataMatrix_Genotype,'wtrescue'),i),'FaceColor',c(4,:),'NumBins',10);
%     h(3,1)=histogram(DataMatrix(strcmp(DataMatrix_Genotype,'het'),i),'FaceColor',c(1,:),'NumBins',10);
%     h(4,1)=histogram(DataMatrix(strcmp(DataMatrix_Genotype,'hetrescue'),i),'FaceColor',c(2,:),'NumBins',10);
%     legend(h,{'wt';'wtrescue';'het';'hetrescue'})

%     h(1,1)=histogram(DataMatrix(strcmp(DataMatrix_Genotype,'wt'),i),'FaceColor',c(3,:),'NumBins',10);
%     h(2,1)=histogram(DataMatrix(strcmp(DataMatrix_Genotype,'het'),i),'FaceColor',c(1,:),'NumBins',10);
%     h(3,1)=histogram(DataMatrix(strcmp(DataMatrix_Genotype,'hetrescue'),i),'FaceColor',c(2,:),'NumBins',10);
      
boxplot([DataMatrix(strcmp(DataMatrix_Genotype,'wt'),i);...
              DataMatrix(strcmp(DataMatrix_Genotype,'hetrescue'),i);...
              DataMatrix(strcmp(DataMatrix_Genotype,'het'),i)],...
              [repmat({['WT (n=',num2str(sum(strcmp(DataMatrix_Genotype,'wt'))),')']},[sum(strcmp(DataMatrix_Genotype,'wt')),1]);...
               repmat({['het+ (n=',num2str(sum(strcmp(DataMatrix_Genotype,'hetrescue'))),')']},[sum(strcmp(DataMatrix_Genotype,'hetrescue')),1]);...
               repmat({['het- (n=',num2str(sum(strcmp(DataMatrix_Genotype,'het'))),')']},[sum(strcmp(DataMatrix_Genotype,'het')),1])],...
               'Colors',c([3,2,1],:))

%       boxplot([DataMatrix(strcmp(DataMatrix_Genotype,'wt'),i);...
%                DataMatrix(strcmp(DataMatrix_Genotype,'wtrescue'),i);...
%                DataMatrix(strcmp(DataMatrix_Genotype,'hetrescue'),i);...
%                DataMatrix(strcmp(DataMatrix_Genotype,'het'),i)],...
%               [repmat({['WT (n=',num2str(sum(strcmp(DataMatrix_Genotype,'wt'))),')']},[sum(strcmp(DataMatrix_Genotype,'wt')),1]);...
%                repmat({['WT+ (n=',num2str(sum(strcmp(DataMatrix_Genotype,'wtrescue'))),')']},[sum(strcmp(DataMatrix_Genotype,'wtrescue')),1]);...
%                repmat({['het+ (n=',num2str(sum(strcmp(DataMatrix_Genotype,'hetrescue'))),')']},[sum(strcmp(DataMatrix_Genotype,'hetrescue')),1]);...
%                repmat({['het- (n=',num2str(sum(strcmp(DataMatrix_Genotype,'het'))),')']},[sum(strcmp(DataMatrix_Genotype,'het')),1])],...
%                'Colors',c([3,4,2,1],:))
           
%     h(1,1)=histogram(DataMatrix(strcmp(DataMatrix_Genotype,'wt'),i),'FaceColor',c(3,:));
%     h(2,1)=histogram(DataMatrix(strcmp(DataMatrix_Genotype,'het'),i),'FaceColor',c(1,:),'BinWidth',get(h(1),'BinWidth'));
%     legend(h,{'wt';'het'})
    
    xlabel('Genotype','Interpreter','latex')
    
    switch i
        case 1
            ylabel('Percentage Awake of total Light period (\%)','Interpreter','latex')
            title('Awake Light percentage','Interpreter','latex')
        case 2
            ylabel('Percentage NREM of total Light period (\%)','Interpreter','latex')
            title('NREM Light percentage','Interpreter','latex')
        case 3
            ylabel('Percentage REM of Sleep in Light period (\%)','Interpreter','latex')
            title('REM Sleep in Light percentage','Interpreter','latex')
        case 4
            ylabel('Percentage Awake of total Dark period (\%)','Interpreter','latex')
            title('Awake Dark percentage','Interpreter','latex')
        case 5
            ylabel('Percentage NREM of total Dark period (\%)','Interpreter','latex')
            title('NREM Dark percentage','Interpreter','latex')
        case 6
            ylabel('Percentage REM of Sleep in Dark period (\%)','Interpreter','latex') 
            title('REM Sleep in Dark percentage','Interpreter','latex')
        case 7
            ylabel('Mean time length (s)','Interpreter','latex')
            title('Means of Awake-light time length distributions for different mice','Interpreter','latex')
        case 8
            ylabel('Mean time length (s)','Interpreter','latex')
            title('Means of NREM-light time length distributions for different mice','Interpreter','latex')
        case 9
            ylabel('Mean time length (s)','Interpreter','latex')
            title('Means of REM-light time length distributions for different mice','Interpreter','latex')
        case 10
            ylabel('Mean time length (s)','Interpreter','latex')
            title('Means of Awake-dark time length distributions for different mice','Interpreter','latex')
        case 11
            ylabel('Mean time length (s)','Interpreter','latex')
            title('Means of NREM-dark time length distributions for different mice','Interpreter','latex')
        case 12
            ylabel('Mean time length (s)','Interpreter','latex')
            title('Means of REM-dark time length distributions for different mice','Interpreter','latex')
        case 13
            ylabel('STD time length (s)','Interpreter','latex')
            title('STD of Awake-light time length distributions for different mice','Interpreter','latex')
        case 14
            ylabel('STD time length (s)','Interpreter','latex')
            title('STD of NREM-light time length distributions for different mice','Interpreter','latex')
        case 15
            ylabel('STD time length (s)','Interpreter','latex')
            title('STD of REM-light time length distributions for different mice','Interpreter','latex')
        case 16
            ylabel('STD time length (s)','Interpreter','latex')
            title('STD of Awake-dark time length distributions for different mice','Interpreter','latex')
        case 17
            ylabel('STD time length (s)','Interpreter','latex')
            title('STD of NREM-dark time length distributions for different mice','Interpreter','latex')
        case 18
            ylabel('STD time length (s)','Interpreter','latex')
            title('STD of REM-dark time length distributions for different mice','Interpreter','latex')
    end
    
    Data_Current={DataMatrix(strcmp(DataMatrix_Genotype,'wt'),i);...
                  DataMatrix(strcmp(DataMatrix_Genotype,'het'),i);...
                  DataMatrix(strcmp(DataMatrix_Genotype,'hetrescue'),i)};
              
    for g=1:3
    for gg=1:3
        [H{i}(g,gg),P{i}(g,gg),CI{i}{g,gg},STATS{i}{g,gg}]=ttest2(Data_Current{g},Data_Current{gg},'vartype','unequal');
    end
    end
    P_Significance(i)=round(sum(sum(P{i}<p_threshold))/2);
    
end
P_Significance=[P_Significance,(1:18)'];

%% NB:
% suggestion needed to be removed for '15332-05_het_D1_Analyzed.mat'
%      15332-03 3/4 right-M/right-S CC=.9962

%% Convert DataTable5 to .xlsx
% DataTable_Current=DataTable5;
% [{'Name'},{'Genotype'},{'Day'},{'State'},{'Channel'},num2cell(fff);DataTable_Current.Name,DataTable_Current.Genotype,DataTable_Current.Day,DataTable_Current.State,DataTable_Current.Channel,num2cell([DataTable_Current.PSD{:}]')]
% xlswrite('AS_PSD.xlsx',ans)
%% Convert DataTable3 to .xlsx
% DataTable_Current=DataTable3;
% for i=1:size(DataTable_Current,1)
% temP(i)=length(DataTable_Current.LengthDistribution{i});
% end
% TEMP=cell(size(DataTable_Current,1),max(temP));
% for i=1:size(DataTable_Current,1)
% TEMP(i,1:temP(i))=num2cell(DataTable_Current.LengthDistribution{i}');
% end
% [{'Name'},{'Genotype'},{'Day'},{'State'},{'LengthDistribution'},cell(1,max(temP));DataTable_Current.Name,DataTable_Current.Genotype,DataTable_Current.Day,DataTable_Current.State,DataTable_Current.LengthDistribution,TEMP]
% xlswrite('AS_LengthDistributions.xlsx',ans)

%% In the case that the 2nd decision failed
%     if sum(logical.Name&strcmp(State_cell,'REM_Light'))==0

%% Save all figures
% FolderName = 'C:\Users\enzo\Downloads';   % Your destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
% FigHandle = FigList(iFig);
% FigName   = get(FigHandle, 'Name');
% saveas(FigHandle, fullfile(FolderName, [FigName,'.png']));
% end
%% Zoom in on delta/theta
% for i=1:6
% subplot(3,2,i)
% xlim([0.5   11.7088])
% end