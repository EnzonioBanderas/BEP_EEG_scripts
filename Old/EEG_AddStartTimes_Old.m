clearvars,close all
EF = uigetdir('','Select Experiment Folder (EF)');
cd(EF)
E=dir('*_Analyzed.mat');
nE=length(E);
OPTIONS.Resize='on';
OPTIONS.WindowStyle='normal';

if exist('DataTable_StartTime.mat','file')
    
    load('DataTable_StartTime')
    count=0;
    for i=1:nE
        if  ~any(strcmp(DataTable_StartTime.Name,E(i).name))
            count=count+1;
            E_Name_inputcell{count,1}=E(i).name;
        end
    end
    
    if count>0
        starttime_default=cell(count,1);
        starttime_default(:)={'09.00.00'};
        starttimes=inputdlg(E_Name_inputcell,'Enter start times (e.g. 09.49.57)',1,starttime_default,OPTIONS);
        DataTable_StartTime=[DataTable_StartTime;table(E_Name_inputcell,starttimes,'VariableNames',{'Name';'StartTime'})];
        save('DataTable_StartTime','DataTable_StartTime')
    else
        msgbox('all experiments already have starting times')
    end
    
else
    
    E_Name_inputcell={E.name}';
    starttime_default=cell(nE,1);
    starttime_default(:)={'09.00.00'};
    starttimes=inputdlg(E_Name_inputcell,'Enter start times (e.g. 09.49.57)',1,starttime_default,OPTIONS);
    DataTable_StartTime=table(E_Name_inputcell,starttimes,'VariableNames',{'Name';'StartTime'});
    save('DataTable_StartTime','DataTable_StartTime')
    
end



