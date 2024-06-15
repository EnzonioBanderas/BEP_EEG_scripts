% clearvars, close all
% Select to be analyzed .mat files
[E_Name,PathName] = uigetfile('*Analyzed.mat','Select the file to analyse','MultiSelect', 'on');
if ischar(E_Name) %if-statement for the case that only one file is selected
    E_Name={E_Name}; 
end
cd(PathName)
nE=length(E_Name);

load('Settings')

figure; colorPlot=colormap(lines); close(gcf)
colorPlot([1,2],:)=colorPlot([2,1],:);

CC_mat=cell(nE,1);
check=triu(true(Settings.nChannel),1);
for i=1:nE
    
    CC_mat{i}=ones(Settings.nChannel);
    
    load(E_Name{i},'Epoch','PSD','FFT')
    
    % compute
    for ROW=1:Settings.nChannel
    for COL=1:Settings.nChannel
    if check(ROW,COL)
        CC_mat{i}(ROW,COL)=mean(Epoch.CC{ROW,COL}(:,end));
    end
    end
    end
    
    figure('name',E_Name{i},'windowstate','maximized')
    
        % left-S
        subplot(2,3,2), hold on, count=0;
        for iii=1:Settings.nState
            if ~isempty(PSD{4,iii})
                count=count+1;
                plot(FFT.BinFreq,PSD{4,iii},'LineWidth',2,'color',colorPlot(iii,:))
            end
        end
        title(Settings.Channels{4})
        xlim([0,30])
        legend(Settings.State(1:count))
        grid on
        % right-S
        subplot(2,3,5), hold on, count=0;
        for iii=1:Settings.nState
            if ~isempty(PSD{1,iii})
                count=count+1;
                plot(FFT.BinFreq,PSD{1,iii},'LineWidth',2,'color',colorPlot(iii,:))
            end
        end
        title(Settings.Channels{1})
        xlim([0,30])
        legend(Settings.State(1:count))
        grid on
        % left-M
        subplot(2,3,1), hold on, count=0;
        for iii=1:Settings.nState
            if ~isempty(PSD{3,iii})
                count=count+1;
                plot(FFT.BinFreq,PSD{3,iii},'LineWidth',2,'color',colorPlot(iii,:))
            end
        end
        title(Settings.Channels{3})
        xlim([0,30])
        legend(Settings.State(1:count))
        grid on
        % right-M
        subplot(2,3,4), hold on, count=0;
        for iii=1:Settings.nState
            if ~isempty(PSD{2,iii})
                count=count+1;
                plot(FFT.BinFreq,PSD{2,iii},'LineWidth',2,'color',colorPlot(iii,:))
            end
        end
        title(Settings.Channels{2})
        xlim([0,30])
        legend(Settings.State(1:count))
        grid on
        
    subplot(2,3,3)
    imagesc(CC_mat{i})
    colorbar
    
end