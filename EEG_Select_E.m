function select_points = EEG_Select_E(data,data_t,data_lin,EF_Name,Settings,Channels)
% LFP_Select_EF takes data processed by LFP_Process_EF and plots a
% selection window where a selection of data can be made by clicking on the
% figure. The first click will mark the start of the first selection
% segment, the second click the end, after which this process can be repeated
% for an arbitrary amount of segments. The time threshold for when a
% selection is accepted is defined in LFP_Process_EF. If the selection time
% is greater than the threshold the title shows up as green, if the
% selection time is lower than the threshold the title shows up as red. The
% selection is converted to a nSegment (number of segments) by 2 matrix, with the 
% first column containing starting points of segments and the second column 
% containing end points of segments.
%
% Enzo Nio, Sep 2018
% enzonio@hotmail.com
%% Select data
nChannel=size(Channels,1);
SelectMultipleChannels=Settings.SelectMultipleChannels;
ShowLinearFit=Settings.ShowLinearFit;
Randomize=Settings.Randomize;

LFP_STDmultiplier=5;

% % Seperate experiment folder name
% if ~Randomize
% EF_Name = strsplit(EF_Name, '_');
% EF_Name = [EF_Name{1}, '_', EF_Name{4}];
% end
    
% % Create new figure if figure does not exist, else select current figure
% g = groot;
% if isempty(g.CurrentFigure) % True if there is no current figure
%     select_fig=figure('units','normalized','outerposition',[0 0 1 1]);
% else
%     select_fig=gcf;
% end

% Create new figure
select_fig=figure('units','normalized','outerposition',[0 0 1 1]);

% Preallocate
select_points=[];

% Colors for plotting
c=[0,0.4470,0.7410;0.8500,0.3250,0.0980;0.4660,0.6740,0.1880]; % c=[blue;red;green]

% % Suggestion for selection
% suggestion=true(length(data_t),1);
% kernel=true(5*Settings.fs,1); % Kernel of 5 s
% for i=1:nChannel
%     suggestion=suggestion&~conv(data{i}>3*std(data{i}),kernel,'same');
% end
% suggestion=logical2points(suggestion);
% nSuggestion=size(suggestion,1);
    
while isempty(select_points)

select_points=cell(nChannel,1);
    
% Plot channels in one figure if input for Settings.SelectMultipleChannels is true
if SelectMultipleChannels
    
    % Name figure with experiment folder name input
    select_fig.Name=['Select time window for ',EF_Name];
    
    % Plot all channels with suggested selection in one figure 
    for i = 1:nChannel
        subplot(nChannel,1,i)
        hold on
        grid on
        plot(data_t,data(:,i),'Color',c(1,:))
%         for ii=1:nSuggestion
%             index=suggestion(ii,1):suggestion(ii,2);
%             plot(data_t(index),data{i}(index),'Color',c(3,:))
%         end
        % Plot linear fit to unfiltered data if Settings.ShowLinearFit is 'on'
        if ShowLinearFit
            plot(data_t,data_lin{i},'Color',c(2,:))
        end
        xlim([0,data_t(end)])
%         YLIM(1)=-LFP_STDmultiplier*std(data(:,i));
%         YLIM(2)=-YLIM(1);
%         ylim(YLIM)
        xlabel('Time (s)')
        if i==5
            ylabel('Acceleration (g)')   
        else
            ylabel('EEG (\muV)')   
        end
%         suggested_time=sum(suggestion(:,2)-suggestion(:,1)+1)/Settings.fs;
%         if suggested_time>Settings.time_threshold*60
%             title([Channels{i,1},' ',num2str(suggested_time),' s'],'Color',c(3,:),'Interpreter','none')
%         else
%             title([Channels{i,1},' ',num2str(suggested_time),' s'],'Color',c(2,:),'Interpreter','none')
%         end
        title(Channels{i,1})
    end
        
    % Select time frame
    m=msgbox('Select time frames?');
    waitfor(m)
    pointx_temp=0; index=0; 
%     selected_time=0;
    while ~isempty(pointx_temp)
        
        [pointx_temp,~,button]=ginput(1);
        
        % If selected point is to the left or right of the axis it will be 
        % assigned to the start or end of the LFP trace respectively.
        if pointx_temp<data_t(1)
            pointx_temp=data_t(1);
        elseif pointx_temp>data_t(end)
            pointx_temp=data_t(end);
        end
        
        % If 'p' key is pressed during selection, wait for the pressing of the
        % 'c' key and then continue selecting.
        if button==double('p')
            waitfor(gcf,'CurrentCharacter','c')
        
        % If 'r' key is pressed during selection the selection is restarted
        % by removing previously selected points
        elseif button==double('r')
            index=0;
            select_points_temp=[];
        
        % If any other key than 'p' or 'Enter' is pressed, the x coordinate
        % of the ginput output will not be empty and will be used to alter
        % the figure and be saved in the selection.
        elseif ~isempty(pointx_temp)
    
            index=index+1;
            % plot over suggestion after first click
            if index==1
                for i=1:nChannel
                    subplot(nChannel,1,i)
                    cla
                    plot(data_t,data(:,i),'Color',c(1,:))
                end
            end
            
            select_points_temp(index,1)=round(pointx_temp*Settings.fs)+1;
            for i=1:nChannel % plot for all channels
                subplot(nChannel,1,i)
                plot(data_t(select_points_temp(index)),data(select_points_temp(index),i),'rx')
            end
            
            if mod(index,2)==0
                for i=1:nChannel %plot for all channels and use same points
                % plot
                subplot(nChannel,1,i)
                plot_index=select_points_temp(index-1):select_points_temp(index);
                plot(data_t(plot_index),data(plot_index,i),'Color',c(3,:))
                
                
                selected_time=sum(points2logical([select_points_temp(1:2:end),select_points_temp(2:2:end)]))/Settings.fs;
                if selected_time>Settings.time_threshold*60
                    title([Channels{i,1},' ',num2str(selected_time),' s'],'Color',c(3,:),'Interpreter','none')
                else
                    title([Channels{i,1},' ',num2str(selected_time),' s'],'Color',c(2,:),'Interpreter','none')
                end
                
                end
            end
            
        end
        
    end % while loop ends when 'Enter' key is pressed
            
    if index==0 % no points selected
        for i=1:nChannel
            select_points{i}=[];
        end
    else
        nSegment=floor(length(select_points_temp)/2);
        select_points_temp=[select_points_temp(1:2:nSegment*2),select_points_temp(2:2:nSegment*2)];
        select_logical_temp=points2logical(select_points_temp);
        for i=1:nChannel        
            select_points{i}=logical2points(select_logical_temp);
        end
    end
        
    clf
    
else % Settings.SelectMultipleChannels is false
    
    for i=1:nChannel
        
        % Reset select_points_temp after first channel
        if i>1
            select_points_temp=[];
        end
        
        % Name figure with experiment folder name and channel name inputs
        select_fig.Name=['Select time window for ',EF_Name,' ',Channels{i,1}];
    
        hold on
        grid on
        plot(data_t,data(:,i),'Color',c(1,:))
%         % Plot suggestion
%         for ii=1:nSuggestion
%             index=suggestion(ii,1):suggestion(ii,2);
%             plot(data_t(index),data{i}(index),'Color',c(3,:))
%         end
        % Plot linear fit to unfiltered data if Settings.ShowLinearFit is 'on'
        if ShowLinearFit
            plot(data_t,data_lin{i},'Color',c(2,:))
        end
        xlim([0,data_t(end)])
%         YLIM(1)=-LFP_STDmultiplier*std(data(:,i));
%         YLIM(2)=-YLIM(1);
%         ylim(YLIM)
        xlabel('Time (s)')
        if i==5
            ylabel('Acceleration (g)')   
        else
            ylabel('EEG (\muV)')   
        end
%         suggested_time=sum(suggestion(:,2)-suggestion(:,1)+1)/Settings.fs;
%         if suggested_time>Settings.time_threshold*60
%             title([Channels{i,1},' ',num2str(suggested_time),' s'],'Color',c(3,:),'Interpreter','none')
%         else
%             title([Channels{i,1},' ',num2str(suggested_time),' s'],'Color',c(2,:),'Interpreter','none')
%         end
        title(Channels{i,1})
        
        % Select time frame
        m=msgbox('Select time frames?');
        waitfor(m)
        pointx_temp=0; index=0; 
%         selected_time=0;
        while ~isempty(pointx_temp)
            
            [pointx_temp,~,button]=ginput(1);
        
        % If selected point is to the left or right of the axis it will be 
        % assigned to the start or end of the LFP trace respectively.
        if pointx_temp<data_t(1)
            pointx_temp=data_t(1);
        elseif pointx_temp>data_t(end)
            pointx_temp=data_t(end);
        end
        
        % If 'p' key is pressed during selection, wait for the pressing of the
        % 'c' key and then continue selecting.
        if button==double('p')
            waitfor(gcf,'CurrentCharacter','c')
            
        % If 'r' key is pressed during selection the selection is restarted
        % by removing previously selected points
        elseif button==double('r')
            index=0;
            select_points_temp=[];
        
        % If any other key than 'p' or 'Enter' is pressed, the x coordinate
        % of the ginput output will not be empty and will be used to alter
        % the figure and be saved in the selection.
        elseif ~isempty(pointx_temp)
    
                index=index+1;
                if index==1
                    plot(data_t,data(:,i),'Color',c(1,:))
                    title(Channels{i,1},'Color',c(2,:),'Interpreter','none')
                end
                select_points_temp(index,1)=round(pointx_temp*Settings.fs)+1;
                plot(select_points_temp(index),data(select_points_temp(index),i),'rx')
            
                if mod(index,2)==0
                    plot_index=select_points_temp(index-1):select_points_temp(index);
                    % plot
                    plot(data_t(plot_index),data(plot_index,i),'Color',c(3,:))
                
                selected_time=sum(points2logical([select_points_temp(1:2:end),select_points_temp(2:2:end)]))/Settings.fs;
                if selected_time>Settings.time_threshold*60
                    title([Channels{i,1},' ',num2str(selected_time),' s'],'Color',c(3,:),'Interpreter','none')
                else
                    title([Channels{i,1},' ',num2str(selected_time),' s'],'Color',c(2,:),'Interpreter','none')
                end
                
                end
                
        end
        end % while loop ends when 'Enter' key is pressed
            
        if index==0 % no points selected
            select_points{i}=[];
        else
            nSegment=floor(length(select_points_temp)/2);
            select_points_temp=[select_points_temp(1:2:nSegment*2),select_points_temp(2:2:nSegment*2)];
            select_logical_temp=points2logical(select_points_temp);
            select_points{i}=logical2points(select_logical_temp);
        end
        
        clf
    
    end
    
end

%% Adjust select_points to only include segments which have a point length greater than the window length
% >=nWindow
for i=1:nChannel
    select_points{i}=select_points{i}(select_points{i}(:,2)-select_points{i}(:,1)+1>=floor(Settings.fs*Settings.window),:);
end

%% If there are no selected segments left for at least one of the channels restart selection
select_points_check=false(nChannel,1);
for i=1:nChannel
    select_points_check(i)=isempty(select_points{i});
end
if any(select_points_check)
    select_points=[];
    msgbox('Not enough points selected')
end

end


end