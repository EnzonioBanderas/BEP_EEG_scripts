clearvars,close all
EF = uigetdir('','Select Experiment Folder (EF)');
cd(EF)
E=dir('*_Analyzed.mat');
nE=length(E);
OPTIONS.Resize='on';
OPTIONS.WindowStyle='normal';

% if DataTable containing starting times of measurements already exists,
% then load DataTable and check whether any of the experiments in the
% selected folder already have starting times.
if exist('DataTable_StartTime.mat','file')
    
    load('DataTable_StartTime')
    count=0;
    for i=1:nE
        if  ~any(strcmp(DataTable_StartTime.Name,E(i).name))
            count=count+1;
            E_Name_inputcell{count,1}=E(i).name;
        end
    end
    
    if count>0 % if at least one of the experiments does not have a starting time
    starttime_default=cell(count,1);
    starttime_default(:)={'2018-10-15 09:00:00'};
    starttimes=inputdlg(E_Name_inputcell,'Enter start times (e.g. 2017-09-16 09:49:57)',1,starttime_default,OPTIONS);
    % datetime of each state
    for i=1:count
        starttimes{i}=datetime(starttimes{i});
    end
    if count==1
        DataTable_StartTime=[DataTable_StartTime;table(E_Name_inputcell,starttimes{:},'VariableNames',{'Name';'StartTime'})];
    else
        DataTable_StartTime=[DataTable_StartTime;table(E_Name_inputcell,cat(1,starttimes{:}),'VariableNames',{'Name';'StartTime'})];
    end
    save('DataTable_StartTime','DataTable_StartTime')
    else % if all experiments already have starting times
        msgbox('all experiments already have starting times')
    end
    
else % if no DataTable exists within selected folder
    
    E_Name_inputcell={E.name}';
    starttime_default=cell(nE,1);
    starttime_default(:)={'2018-10-15 09:00:00'};
    starttimes=inputdlg(E_Name_inputcell,'Enter start times (e.g. 2017-09-16 09:49:57)',1,starttime_default,OPTIONS);
    % datetime of each state
    for i=1:nE
        starttimes{i}=datetime(starttimes{i});
    end
    if nE==1
        DataTable_StartTime=table(E_Name_inputcell,starttimes{:},'VariableNames',{'Name';'StartTime'});
    else
        DataTable_StartTime=table(E_Name_inputcell,cat(1,starttimes{:}),'VariableNames',{'Name';'StartTime'});
    end
    save('DataTable_StartTime','DataTable_StartTime')
    
end



