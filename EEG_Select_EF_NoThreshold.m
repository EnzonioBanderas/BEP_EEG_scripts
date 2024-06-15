% Select part of trace that is to be used in analysis.
% For the selection of a different part the old DataTable_Select.mat
% should first be removed
clearvars, close all
% select directory
experiment_folder = uigetdir('','Select folder containing experiments');
cd (experiment_folder)
experiments = dir('*_*_*.mat');
nExperiments = length(experiments);

% Settings
Settings.timeBin=60; % s
COLOR  =   [0         0.4470    0.7410
            0.8500    0.3250    0.0980];
YLIM=[-1500,1500];

count=0;
for i=1:nExperiments

load(experiments(i).name)
nChannel=size(Electrical.CH1234,2);
Acceleration.t=0:1/Acceleration.fs:(size(Acceleration.XYZ,1)-1)/Acceleration.fs;
Electrical.t=0:1/Electrical.fs:(size(Electrical.CH1234,1)-1)/Electrical.fs;
    
check=true;
while check
    
    % if experiment select is not already there
    count=count+1;
        
    load(experiments(i).name,'Acceleration','Electrical')
    nChannels=size(Electrical.CH1234,2);
    Acceleration.t=0:1/Acceleration.fs:(size(Acceleration.XYZ,1)-1)/Acceleration.fs;
    Electrical.t=0:1/Electrical.fs:(size(Electrical.CH1234,1)-1)/Electrical.fs;
    
    % Mouse code
    Mouse = experiments(i).name(1:end-4);
    Mouse(experiments(i).name(1:end-4)=='_')=' ';

    % Plot data       
    figure('Name',['Select time window for ', Mouse],...
           'units','normalized','outerposition',[0 0 1 1])
    grid on
    for ii=1:nChannels
        subplot(nChannels,1,ii)
        plot(Electrical.t,Electrical.CH1234(:,ii))
        xlim([0,Electrical.t(end)]), ylim([-1000,1000])
        xlabel('Time (s)'), ylabel('Field Potential (\muV)')
    end
        
    % Select time frame
    m=msgbox('Select time frame?');
    waitfor(m)
    select_logical_acc=false(1,length(Acceleration.t));
    select_logical_ele=false(1,length(Electrical.t));
    pointx_temp=0; index=0;
    while ~isempty(pointx_temp)
        index=index+1;
        [pointx_temp,~]=ginput(1);
        if ~isempty(pointx_temp)
            pointx(index)=pointx_temp;
            for ii=1:nChannels
                subplot(nChannels,1,ii), hold on
                plot(ones(2,1)*pointx(index),YLIM,'r','LineWidth',2)
            end
            if mod(index,2)==0
                % Electrical Logical
                select_logical_temp=Electrical.t>=pointx(index-1)&...
                                    Electrical.t<=pointx(index);
                select_logical_ele=select_logical_ele|select_logical_temp;
                % Acceleration logical
                select_logical_temp=Acceleration.t>=pointx(index-1)&...
                                    Acceleration.t<=pointx(index);
                select_logical_acc=select_logical_acc|select_logical_temp;
            end 
        end
    end
    close;
        
    experiments_cell{count,1}=experiments(i).name(1:end-4);
    select_points_acc_cell{count,1}=logical2points(select_logical_acc);
    select_points_ele_cell{count,1}=logical2points(select_logical_ele);
    
    select_points_acc=logical2points(select_logical_acc);
    select_points_ele=logical2points(select_logical_ele);
    
    % remove selections that are less than bin length
    select_points_acc=select_points_acc(select_points_acc(:,2)-select_points_acc(:,1)+1>=...
    Settings.timeBin*Acceleration.fs,:);
    select_points_ele=select_points_ele(select_points_ele(:,2)-select_points_ele(:,1)+1>=...
    Settings.timeBin*Electrical.fs,:);

    if ~(isempty(select_points_acc)||isempty(select_points_ele))
        check=false;
    end
    
end

Acceleration.points=select_points_acc;
EleAccMultiple=round(Electrical.fs/Acceleration.fs);
Electrical.points=Acceleration.points*EleAccMultiple;
Electrical.points(:,1)=Electrical.points(:,1)-EleAccMultiple+1;
save(experiments(i).name,'Acceleration','Electrical','-append')

end