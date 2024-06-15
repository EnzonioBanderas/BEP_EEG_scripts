% Plot mean PSD for each mouse and genotype
clearvars, close all
%% Import
EF = uigetdir('','Select Experiment Folder (EF)');
cd(EF)
E=dir('*_Analyzed.mat');
nE=length(E);
load('DataTable')

load('Settings')

% Leave out high cross correlation between channels
LeaveOut=strcmp(DataTable4.Name,'14248-05')&(strcmp(DataTable4.Channel,'right-M')|strcmp(DataTable4.Channel,'right-S'))|...
        strcmp(DataTable4.Name,'14346-01')&(strcmp(DataTable4.Channel,'left-S') |strcmp(DataTable4.Channel,'left-M')) |...
        strcmp(DataTable4.Name,'14346-03')&(strcmp(DataTable4.Channel,'right-M')|strcmp(DataTable4.Channel,'right-S')|strcmp(DataTable4.Channel,'right-S') |strcmp(DataTable4.Channel,'right-M'))|...
        strcmp(DataTable4.Name,'14351-02')&(strcmp(DataTable4.Channel,'right-S')|strcmp(DataTable4.Channel,'right-M'))|...
        strcmp(DataTable4.Name,'14463-02')&(strcmp(DataTable4.Channel,'right-M')|strcmp(DataTable4.Channel,'right-S'));
DataTable4(LeaveOut,:)=[];

%% Constants defined here (can be changed to input prompt)
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


uniq.Name_full=unique(DataTable4.Name);
uniq.State_full=unique(DataTable4.State);
nState_full=length(uniq.State_full);

count=0;
uniq.Channel=unique(DataTable4.Channel);
nChannel=length(uniq.Channel);
for i=1:nChannel
    logical.channel=strcmp(DataTable4.Channel,uniq.Channel(i));
    
    figure('Name',[uniq.Channel{i},' Mean over Days'])
    legend_cell=cell(nState_full,2);
    handle_cell=cell(nState_full,2);
    
    uniq.Name=unique(DataTable4.Name(logical.channel));
    
    for ii=1:length(uniq.Name)
    logical.channel_name=logical.channel&strcmp(DataTable4.Name,uniq.Name{ii});
    uniq.State=unique(DataTable4.State(logical.channel_name));
    
    for iii=1:length(uniq.State)
    logical.channel_name_state=logical.channel_name&strcmp(DataTable4.State,uniq.State{iii});
        
    if sum(logical.channel_name_state)
    
    count=count+1;    
        
    % do mean over days and save in cell for temp table
    PSD_Matrix=[DataTable4.PSD{logical.channel_name_state}];
    [temptable.PSD{count,1},SEM_meanday]=meanSEM(PSD_Matrix);
    
    % additional strings for temp table
    temptable.Channel(count,1)=uniq.Channel(i);
    temptable.Genotype(count,1)=DataTable4.Genotype(find(logical.channel_name_state,1,'first'));
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
    subplot(nState_full,2,(state_index-1)*2+1) % index for plot (2x2 : state x norm)
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
    legend_cell{iii,1}=[legend_cell{iii,1},DataTable4.Name(find(logical.channel_name_state,1,'first'))];
    handle_cell{iii,1}=[handle_cell{iii,1},h];
    xlabel('Frequency (Hz)'), ylabel('PSD (\muV^2)')
    
    % normalization (%)
    PSD_Matrix=(PSD_Matrix./bandpower2(PSD_Matrix,fff,[Settings.HP_ele,normFreq]))*100; 
    [PSD_meanday_normalized,SEM_meanday_normalized]=meanSEM(PSD_Matrix);
    
    % normalized plot
    subplot(nState_full,2,state_index*2) % index for plot
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
    legend_cell{iii,2}=[legend_cell{iii,2},DataTable4.Name(find(logical.channel_name_state,1,'first'))];
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
        LEG=legend(handle_cell{1,1}',legend_cell{1,1}');
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
    nState_full=length(uniq.State);
    for ii=1:nState_full
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
    

% for l=1:nState_full
%     for ll=1:2
%         subplot(nState_full,2,ll+(l-1)*2)
%         legend(handle_cell{l,ll},legend_cell{l,ll})
%     end
% end

        subplot(nState_full,2,1)
        LEG=legend(handle_cell{1},legend_cell{1});
        LEG.Position=[0.4807 0.4929 0.0608 0.0578];
    
end 