clearvars
close all

% EF = uigetdir('','Select Experiment Folder (EF)');
% cd(EF)
% E=dir('*_Analyzed.mat');
% nE=length(E);

% Select to be analyzed .mat files
[E_Name,PathName] = uigetfile('*Analyzed.mat','Select the file to analyse','MultiSelect', 'on');
if ischar(E_Name) %if-statement for the case that only one file is selected
    E_Name={E_Name}; 
end
cd(PathName)
nE=length(E_Name);
for i=1:nE
    E(i).name=E_Name{i};
end

load('DataTable_StartTime')

% Interpolation: When the state is unknown interpolate to surrounding state
% if the length of unknown segment is also under appropriate defined interpolation time.
AwakeInterpolationTime=3*60;
 NREMInterpolationTime=3*60;
  REMInterpolationTime=30;

% Rule 1: QA between REM or NREM and REM is likely also REM (<REMInterpolationTime)
REM2QAInterpolationTime=30;
% Rule 2: REM does not follow Awake(=AA||QA), instead this is probably QA. Only AA should be used (QA untrustworthy). 
REMafterAwakeCutoffTime=3*60;

% Get sampling rate and assume to it be constant across measurements
load(E(1).name,'Electrical')
fs=Electrical.fs;
clearvars Electrical

c=colormap(lines);
c([2,1],:)=c([1,2],:);
close(gcf)

% Segment length for circadian statistics
segmentLength=round(60*60*fs);
segmentOverlap=round(segmentLength*(9/10));
dt=(segmentLength-segmentOverlap)/fs;

nowtrescue=false;
noalign=true;

% Define between which hours it is light and dark
Light=[7.5,18.5];
Dark=[19.5,6.5];
% maxTime=96; % hours
stateNames={'Awake';'NREM';'REM'};
nstateNames=length(stateNames);
genotypeNames={'het';'wt';'hetrescue';'wtrescue'};
if nowtrescue
    genotypeNames(end)=[];
end
ngenotypeNames=length(genotypeNames);
% cut_off=0.5*10^-4;
cut_off=2.5*10^-5;

% nE=2;
circadianTime=cell(nE,1);
circadianStatistics=cell(nE,1);
StatePointsCell=cell(nE,1);
circadianFiltCell=cell(nE,1);
circadianFiltTime=cell(nE,1);
fs_circadian=zeros(nE,1);
nElectrical=zeros(nE,1);
rejectE=false(nE,1);
nInterp=zeros(nE,1);
nQAInterp=zeros(nE,1);
nREM2QA=zeros(nE,1);
for i=1:nE
    
    load(E(i).name,'Electrical_P')
    nanLogical=~isnan(Electrical_P(:,1));
    [~,state_index]=max(Electrical_P(nanLogical,:),[],2);
    clearvars Electrical_P
    nanLogical_double=double(nanLogical);
    nanLogical_double(nanLogical)=state_index;
    state_index=nanLogical_double;
    clearvars nanLogical_double nanLogical
    
    nElectrical(i)=size(state_index,1);
    
    Unknown=logical2points(state_index==0);
    AA=logical2points(state_index==1);
    NREM=logical2points(state_index==2);
    REM=logical2points(state_index==3);
    QA=logical2points(state_index==4);
    clearvars state_index
    
    nUnknown=size(Unknown,1); % Interpolation
    nQA=size(QA,1); % Rule 1
    nREM=size(REM,1); % Rule 2
    nAA=size(AA,1);
    nNREM=size(NREM,1);
    
    StateID=[repmat({'Unknown'},[nUnknown,1]);...
             repmat({'AA'},[nAA,1]);...
             repmat({'NREM'},[nNREM,1]);...
             repmat({'REM'},[nREM,1]);...
             repmat({'QA'},[nQA,1])];
         
    StatePoints=table([Unknown;AA;NREM;REM;QA],StateID,'variablenames',{'Points','State'});
    StatePoints=sortrows(StatePoints,'Points','ascend');
    
%     % test table, interpolation and rule
%     StatePoints=table([1,100;101,200;201,300;301,400;401,500;501,600;601,700;701,800;801,900;901,1000;1001,1100;1101,1200;1201,1300;1301,1400;1401,1500;1501,1600;1601,1700;1701,1800;1801,1900;1901,2000;2001,2100;2101,2200;2201,2300;2301,2400;2401,2500;2501,2600;2601,2700;2701,2800;2801,2900;2901,3000;3001,3100;3101,3200;3201,3300],...
%     {'Unknown','AA','Unknown','AA','QA','Unknown','QA','NREM','Unknown','NREM','REM','Unknown','REM','QA','REM','AA','REM','AA','Unknown','NREM','Unknown','QA','REM','Unknown','AA','Unknown','REM','NREM','Unknown','REM','Unknown','NREM','Unknown'}','variablenames',{'Points','State'});
%     StatePoints=table([1,1;2,2;3,round(AwakeInterpolationTime*fs)+3;round(AwakeInterpolationTime*fs)+4,round(AwakeInterpolationTime*fs)+4;round(AwakeInterpolationTime*fs)+5,round(AwakeInterpolationTime*fs)+5],...
%                 {'Unknown','AA','Unknown',                              'AA',                                                                   'Unknown'}','variablenames',{'Points','State'});
%     
%     nUnknown=sum(strcmp(StatePoints.State,'Unknown')); % Interpolation
%     nQA=sum(strcmp(StatePoints.State,'QA')); % Rule 1
%     nREM=sum(strcmp(StatePoints.State,'REM')); % Rule 2
    
    nSegment=size(StatePoints.Points,1);
            
    % Interpolation, maybe split up points in halfs and assign state
    UnknownIndex=find(strcmp(StatePoints.State,'Unknown'));
%     nInterp=0;
    for ii=1:nUnknown
        if UnknownIndex(ii)~=1&&UnknownIndex(ii)~=nSegment
            % AA/QA-AA/QA (1,1)
            unknownLength=(StatePoints.Points(UnknownIndex(ii),2)-StatePoints.Points(UnknownIndex(ii),1)+1);
            if (strcmp(StatePoints.State{UnknownIndex(ii)-1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)-1},'QA'))&&...
               (strcmp(StatePoints.State{UnknownIndex(ii)+1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)+1},'QA'))
                if unknownLength/fs<AwakeInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'AA'}; % ignoring QA for now
                    nInterp(i)=nInterp(i)+1;
                end
            % NREM-NREM (2,2)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'NREM')&&...
                   strcmp(StatePoints.State{UnknownIndex(ii)+1},'NREM')
                if unknownLength/fs<NREMInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'NREM'};
                    nInterp(i)=nInterp(i)+1;
                end
            % REM-REM (3,3)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'REM')&&...
                   strcmp(StatePoints.State{UnknownIndex(ii)+1},'REM')
                if unknownLength/fs<REMInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'REM'};
                    nInterp(i)=nInterp(i)+1;
                end
            % AA/QA-NREM (1,2)
            elseif (strcmp(StatePoints.State{UnknownIndex(ii)-1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)-1},'QA'))&&...
                    strcmp(StatePoints.State{UnknownIndex(ii)+1},'NREM')
                if unknownLength/fs<(AwakeInterpolationTime+NREMInterpolationTime)/4
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp(i)=nInterp(i)+1;
                end
            % NREM-AA/QA (2,1)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'NREM')&&...
                  (strcmp(StatePoints.State{UnknownIndex(ii)+1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)+1},'QA'))
                if unknownLength/fs<(AwakeInterpolationTime+NREMInterpolationTime)/4
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp(i)=nInterp(i)+1;
                end
            % AA/QA-REM (1,3)
            elseif (strcmp(StatePoints.State{UnknownIndex(ii)-1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)-1},'QA'))&&...
                   strcmp(StatePoints.State{UnknownIndex(ii)+1},'REM')
                if unknownLength/fs<REMInterpolationTime/2
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp(i)=nInterp(i)+1;
                elseif unknownLength/fs<AwakeInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'AA'}; % ignoring QA
                    nInterp(i)=nInterp(i)+1;
                end
            % REM-AA/QA (3,1)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'REM')&&...
                   (strcmp(StatePoints.State{UnknownIndex(ii)+1},'AA')||strcmp(StatePoints.State{UnknownIndex(ii)+1},'QA'))
                if unknownLength/fs<REMInterpolationTime/2
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp(i)=nInterp(i)+1;
                elseif unknownLength/fs<AwakeInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'AA'}; % ignoring QA
                    nInterp(i)=nInterp(i)+1;
                end
            % NREM-REM (2,3)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'NREM')&&...
                   strcmp(StatePoints.State{UnknownIndex(ii)+1},'REM')
                if unknownLength/fs<REMInterpolationTime/2
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp(i)=nInterp(i)+1;
                elseif unknownLength/fs<NREMInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'NREM'}; % ignoring QA
                    nInterp(i)=nInterp(i)+1;
                end
            % NREM-REM (3,2)
            elseif strcmp(StatePoints.State{UnknownIndex(ii)-1},'REM')&&...
                   strcmp(StatePoints.State{UnknownIndex(ii)+1},'NREM')
                if unknownLength/fs<REMInterpolationTime/2
                    StatePoints.Points(UnknownIndex(ii)-1,2)=StatePoints.Points(UnknownIndex(ii)-1,2)+floor(unknownLength/2);
                    StatePoints.Points(UnknownIndex(ii)+1,1)=StatePoints.Points(UnknownIndex(ii)-1,2)+1;
                    StatePoints.State(UnknownIndex(ii))={'Remove'};
                    nInterp(i)=nInterp(i)+1;
                elseif unknownLength/fs<NREMInterpolationTime
                    StatePoints.State(UnknownIndex(ii))={'NREM'}; % ignoring QA
                    nInterp(i)=nInterp(i)+1;
                end
            end
        end
    end
    StatePoints(strcmp(StatePoints.State,'Remove'),:)=[];
    nSegment=size(StatePoints.Points,1);
    nAA=sum(strcmp(StatePoints.State,'AA'));
    nNREM=sum(strcmp(StatePoints.State,'NREM'));
    nREM=sum(strcmp(StatePoints.State,'REM'));
    nQA=sum(strcmp(StatePoints.State,'QA'));
    nUnknown=sum(strcmp(StatePoints.State,'Unknown'));
    
    % Rule 1: QA between REM or NREM and REM is likely also REM (<REMInterpolationTime)
    QAIndex=find(strcmp(StatePoints.State,'QA'));
%     nQAInterp=0;
    for ii=1:nQA
        if (strcmp(StatePoints.State{QAIndex(ii)-1},'REM')||strcmp(StatePoints.State{QAIndex(ii)-1},'NREM'))&&...
            strcmp(StatePoints.State{QAIndex(ii)+1},'REM')
            if (StatePoints.Points(QAIndex(ii),2)-StatePoints.Points(QAIndex(ii),1)+1)/fs<REM2QAInterpolationTime
                StatePoints.State(QAIndex(ii))={'REM'};
                nQAInterp(i)=nQAInterp(i)+1;
            end
        end
    end
    
    % Rule 2: REM does not follow Awake(=AA||QA), instead this is probably QA. Only AA should be used (QA untrustworthy). 
    REMIndex=find(strcmp(StatePoints.State,'REM'));
%     nREM2QA=0;
%     figure
    for ii=1:nREM
        if REMIndex(ii)~=1
            if strcmp(StatePoints.State{REMIndex(ii)-1},'AA')
                diffREM=(StatePoints.Points(REMIndex(ii),1)-StatePoints.Points(strcmp(StatePoints.State,'AA'),2)+1)/fs;
                if any((diffREM<REMafterAwakeCutoffTime)&(diffREM>0))
                    StatePoints.State(REMIndex(ii))={'QA'};
                    nREM2QA(i)=nREM2QA(i)+1;
                end
%                 plot(diffREM),hold on
            end
        end
    end
    
    % Merge AA and QA to Awake
    StatePoints.State(strcmp(StatePoints.State,'AA')|strcmp(StatePoints.State,'QA'))={'Awake'};
    
    % Merge segments by converting to logical
    UnknownLog=points2logical(StatePoints.Points(strcmp(StatePoints.State,'Unknown'),:),nElectrical(i));
    AwakeLog  =points2logical(StatePoints.Points(strcmp(StatePoints.State,'Awake'),:),nElectrical(i));
    NREMLog   =points2logical(StatePoints.Points(strcmp(StatePoints.State,'NREM'),:) ,nElectrical(i));
    REMLog    =points2logical(StatePoints.Points(strcmp(StatePoints.State,'REM'),:)  ,nElectrical(i));
    UnknownPoints=logical2points(UnknownLog); clearvars UnknownLog
    AwakePoints  =logical2points(AwakeLog);
    NREMPoints   =logical2points(NREMLog);
    REMPoints    =logical2points(REMLog);
    
    % Adjust and store StatePoints table
    StateID=[repmat({'Unknown'},[size(UnknownPoints,1),1]);...
             repmat({'Awake'},[size(AwakePoints,1),1]);...
             repmat({'NREM'},[size(NREMPoints,1),1]);...
             repmat({'REM'},[size(REMPoints,1),1])];
    StatePoints=table([UnknownPoints;AwakePoints;NREMPoints;REMPoints],StateID,'variablenames',{'Points','State'});
    StatePoints=sortrows(StatePoints,'Points','ascend');
    StatePointsCell{i}=StatePoints;
    
    % Calculate circadian statistics
    KnownPoints=logical2points(~points2logical(StatePoints.Points(strcmp(StatePoints.State,'Unknown'),:),nElectrical(i)));   
    [AwakeLog,midPoint]=buffer_dis(AwakeLog,segmentLength,segmentOverlap,KnownPoints); AwakeLog=sum(AwakeLog,1)';
    [NREMLog,~]        =buffer_dis(NREMLog ,segmentLength,segmentOverlap,KnownPoints); NREMLog=sum(NREMLog,1)';
    [REMLog,~]         =buffer_dis(REMLog  ,segmentLength,segmentOverlap,KnownPoints); REMLog=sum(REMLog,1)';
    StartTime=DataTable_StartTime.StartTime(strcmp(DataTable_StartTime.Name,E(i).name));
    circadianTime{i}=midPoint/fs+hour(StartTime)*3600+minute(StartTime)*60+second(StartTime);
    nmidPoint=length(midPoint);
    circadianStatistics{i}=zeros(nmidPoint,3);
    circadianStatistics{i}=[AwakeLog,NREMLog,REMLog];
    circadianStatistics{i}(:,1:3)=(circadianStatistics{i}(:,1:3)./repmat(sum(circadianStatistics{i}(:,1:3),2),[1,3]))*100;
    circadianStatistics{i}(:,3)=(circadianStatistics{i}(:,3)./sum(circadianStatistics{i}(:,2:3),2))*100;
    circadianStatistics{i}(isnan(circadianStatistics{i}(:,3)),3)=0;
    
    % Filter
    fs_circadian=1/dt;
    if ~isempty(circadianTime{i})
    circadianFiltTime{i}=(circadianTime{i}(1):dt:circadianTime{i}(end))';
    circadianFilt=zeros(length(circadianFiltTime{i}),3);
    for ii=1:nstateNames
        %%%%%%%%%%%%%%%%% Would also interpolate >24 hour periods %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        circadianFilt(:,ii)=interp1(circadianTime{i},circadianStatistics{i}(:,ii),circadianFiltTime{i});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        circadianFilt(:,ii)=LowPassfilter(2, cut_off, fs_circadian, circadianFilt(:,ii));
    end
    circadianFiltCell{i}=circadianFilt;
    end
    
    % reject certain experiments
    rejectE(i)=isempty(circadianTime{i});
    
end
% Some experiments will have no 1 hour long segments and should probably be
% rejected for length distribution or circadian analysis
E(rejectE)=[]; nE=length(E); 
StatePointsCell(rejectE,:)=[]; 
circadianFiltCell(rejectE)=[];
circadianFiltTime(rejectE)=[];
circadianStatistics(rejectE)=[];
circadianTime(rejectE)=[];
nElectrical(rejectE)=[];
clearvars -except rejectE E nE StatePointsCell circadianFiltCell circadianFiltTime circadianStatistics circadianTime nElectrical fs dt cut_off c nstateNames stateNames genotypeNames ngenotypeNames DataTable_StartTime  nREM2QA nQAInterp nInterp segmentLength segmentOverlap Light fs_circadian Dark nowtrescue genotypeNames ngenotypeNames noalign nowtrescue

%% Find time delay between mice circadian rhythms
UpDi=triu(true(nE),1);
LoDi=tril(true(nE),-1);
lagMatrix=zeros(nE,nE,nstateNames);
startTimeShift=zeros(nE,nE,nstateNames);
lagVec=zeros(nE,nE);
CCMatrix=zeros(nE,nE,nstateNames);
for iii=1:nstateNames
    
    % calculate cross correlation
    for i=1:nE
    for ii=1:nE
    if UpDi(i,ii)
     
    nMaxlag=round(12*60*60*fs_circadian)*2;
    nWindow=max([length(circadianFiltTime{i}),length(circadianFiltTime{ii})]);
    if nMaxlag>nWindow
        nMaxlag=nWindow;
    end
    [C,~,L] = corrgram2(circadianFiltCell{i}(:,iii) -mean(circadianFiltCell{i}(:,iii)),...
                        circadianFiltCell{ii}(:,iii)-mean(circadianFiltCell{ii}(:,iii)),...
                        nMaxlag,nWindow,0,fs_circadian);
    C=real(C);

    L=L/3600; % hours
    startTimeShift(i,ii,iii)=(circadianFiltTime{ii}(1)-circadianFiltTime{i}(1))/3600; % start time shift in hours
    [pks,locs,~,~]=findpeaks(C,L,'WidthReference','halfheight');
    locs=locs+startTimeShift(i,ii,iii);
    if ~isempty(pks)
        [~,peak_closestto0_index]=min(abs(locs));
        lagMatrix(i,ii,iii)= locs(peak_closestto0_index); % shift i by this value to be aligned to ii
        CCMatrix(i,ii,iii) = pks (peak_closestto0_index);
    else
        lagMatrix(i,ii,iii)= L(nMaxlag+1);
        CCMatrix(i,ii,iii) = C(nMaxlag+1);
    end
%     lagMatrix(i,ii,iii) =lagMatrix(i,ii,iii); +startTimeShift(i,ii,iii)
    
    end    
    end
    end
    
end

% pick row with lowest absolute delay on average
% align all other waves to this wave (column contains necessary shifts)
for i=1:nE
for ii=1:nE
if UpDi(i,ii)
    if CCMatrix(i,ii,1)>=CCMatrix(i,ii,2)
        lagVec(i,ii)=lagMatrix(i,ii,1);
    else
        lagVec(i,ii)=lagMatrix(i,ii,2);
    end
end
end
end
for i=1:nE
for ii=1:nE
    if LoDi(i,ii)
        lagVec(i,ii)=-lagVec(ii,i);
        for iii=1:nstateNames
        lagMatrix(i,ii,iii)=-lagMatrix(ii,i,iii);
        end
    end
end
end
[~,index]=min(sum(abs(lagVec),2)/(nE-1));
lagVec=lagVec(index,:)'; % adjust starting times with this vector (hours)
lagVec2=lagVec;
if noalign
    lagVec=zeros(nE,1);
end

%% use adjusted starting times and exclusion of border segments to calculate new boxplot statistics
StatePoints2Cell=cell(nE,1);
LengthDistribution=cell(nE,6);
DataMatrix_Genotype=cell(nE,1);
DataMatrix_Name=cell(nE,1);
DataMatrix_Day=cell(nE,1);

for i=1:nE
    
    % Genotype
    Esplit=strsplit(E(i).name,'_');
    DataMatrix_Name(i)=Esplit(1);
    DataMatrix_Genotype(i)=Esplit(2);
    DataMatrix_Day(i)=Esplit(3);
    
    % Create Light and Dark logicals
    StartTime=DataTable_StartTime.StartTime(strcmp(DataTable_StartTime.Name,E(i).name));
    StartTime=hour(StartTime)*60*60+minute(StartTime)*60+second(StartTime)-lagVec(i)*3600;
    t=mod((0:1/fs:(nElectrical(i)-1)/fs)'+StartTime,24*60*60);
    LightLog=t>Light(1)*60*60&t<Light(2)*60*60;
    DarkLog=t>Dark(1)*60*60|t<Dark(2)*60*60;
%     DarkLog=~LightLog;
%     DarkPoints=logical2points(DarkLog);
%     DarkLog(DarkPoints(:))=false;
    
    % Retrieve state logicals
    AwakeLog  =points2logical(StatePointsCell{i}.Points(strcmp(StatePointsCell{i}.State,'Awake'),:)  ,nElectrical(i));
    NREMLog   =points2logical(StatePointsCell{i}.Points(strcmp(StatePointsCell{i}.State,'NREM'),:)   ,nElectrical(i));
    REMLog    =points2logical(StatePointsCell{i}.Points(strcmp(StatePointsCell{i}.State,'REM'),:)    ,nElectrical(i));
    
    % Combine Light and Dark logicals with state logicals and turn into points
    AwakeLightPoints=logical2points(AwakeLog&LightLog);
    AwakeDarkPoints =logical2points(AwakeLog&DarkLog);
    NREMLightPoints =logical2points(NREMLog&LightLog);
    NREMDarkPoints  =logical2points(NREMLog&DarkLog);
    REMLightPoints  =logical2points(REMLog&LightLog);
    REMDarkPoints   =logical2points(REMLog&DarkLog);
    clearvars AwakeLog NREMLog REMLog LightLog DarkLog
    StatePointsKnown=[AwakeLightPoints;AwakeDarkPoints;...
                      NREMLightPoints;NREMDarkPoints;...
                      REMLightPoints;REMDarkPoints];
    UnknownPoints=logical2points(~points2logical(StatePointsKnown,nElectrical(i)));
    
    StateID=[repmat({'Unknown'},   [size(UnknownPoints,   1),1]);...
             repmat({'AwakeLight'},[size(AwakeLightPoints,1),1]);...
             repmat({'AwakeDark'}, [size(AwakeDarkPoints, 1),1]);...
             repmat({'NREMLight'}, [size(NREMLightPoints, 1),1]);...
             repmat({'NREMDark'},  [size(NREMDarkPoints,  1),1]);...
             repmat({'REMLight'},  [size(REMLightPoints,  1),1]);...
             repmat({'REMDark'},   [size(REMDarkPoints,   1),1])];
    StatePoints2=table([UnknownPoints;StatePointsKnown],StateID,'variablenames',{'Points','State'});
    StatePoints2=sortrows(StatePoints2,'Points','ascend');
    StatePoints2Cell{i}=StatePoints2;
    
    % Remove segments bordering Unknown
    StatePoints2(conv(double(strcmp(StatePoints2.State,'Unknown')),[1,1,1],'same')>0,:)=[];
    
    LengthDistribution{i,1}=(StatePoints2.Points(strcmp(StatePoints2.State,'AwakeLight'),2)-StatePoints2.Points(strcmp(StatePoints2.State,'AwakeLight'),1)+1)/fs;
    LengthDistribution{i,2}=(StatePoints2.Points(strcmp(StatePoints2.State,'NREMLight'),2)-StatePoints2.Points(strcmp(StatePoints2.State, 'NREMLight') ,1)+1)/fs;
    LengthDistribution{i,3}=(StatePoints2.Points(strcmp(StatePoints2.State,'REMLight'),2)-StatePoints2.Points(strcmp(StatePoints2.State,  'REMLight')  ,1)+1)/fs;
    LengthDistribution{i,4}=(StatePoints2.Points(strcmp(StatePoints2.State,'AwakeDark'),2)-StatePoints2.Points(strcmp(StatePoints2.State, 'AwakeDark') ,1)+1)/fs;
    LengthDistribution{i,5}=(StatePoints2.Points(strcmp(StatePoints2.State,'NREMDark'),2)-StatePoints2.Points(strcmp(StatePoints2.State,  'NREMDark')  ,1)+1)/fs;
    LengthDistribution{i,6}=(StatePoints2.Points(strcmp(StatePoints2.State,'REMDark'),2)-StatePoints2.Points(strcmp(StatePoints2.State,   'REMDark')   ,1)+1)/fs;
    
end
    
% Merge length distributions for same mice and associate genotype with each
DataMatrix_Name_Uniq=unique(DataMatrix_Name);
nName=length(DataMatrix_Name_Uniq);
LengthDistribution2=cell(nName,6);
DataMatrix_Genotype2=cell(nName,1);
% Turn genotype wtrescue into wt
if nowtrescue
    DataMatrix_Genotype(strcmp(DataMatrix_Genotype,'wtrescue'))={'wt'};
end
for i=1:nName
    nameLog=strcmp(DataMatrix_Name,DataMatrix_Name_Uniq{i});
    DataMatrix_Genotype2(i)=DataMatrix_Genotype(find(nameLog,1,'first'));
    for ii=1:6
        LengthDistribution2(i,ii)={cat(1,LengthDistribution{nameLog,ii})};
    end
end
   
DataMatrix=zeros(nName,18);
for i=1:nName
    
    DataMatrix(i,1)=sum(LengthDistribution2{i,1});
    DataMatrix(i,2)=sum(LengthDistribution2{i,2});
    DataMatrix(i,3)=sum(LengthDistribution2{i,3});
    DataMatrix(i,1:3)=[(DataMatrix(i,1:2)/sum(DataMatrix(i,1:3)))*100,(DataMatrix(i,3)/sum(DataMatrix(i,2:3)))*100];
%     DataMatrix(i,1:3)= (DataMatrix(i,1:3)/sum(DataMatrix(i,1:3)))*100;
    DataMatrix(i,4)=sum(LengthDistribution2{i,4});
    DataMatrix(i,5)=sum(LengthDistribution2{i,5});
    DataMatrix(i,6)=sum(LengthDistribution2{i,6});
    DataMatrix(i,4:6)=[(DataMatrix(i,4:5)/sum(DataMatrix(i,4:6)))*100,(DataMatrix(i,6)/sum(DataMatrix(i,5:6)))*100];
%     DataMatrix(i,4:6)= (DataMatrix(i,4:6)/sum(DataMatrix(i,4:6)))*100;
    
    DataMatrix(i,7)=mean( LengthDistribution2{i,1});
    DataMatrix(i,8)=mean( LengthDistribution2{i,2});
    DataMatrix(i,9)=mean( LengthDistribution2{i,3});
    DataMatrix(i,10)=mean(LengthDistribution2{i,4});
    DataMatrix(i,11)=mean(LengthDistribution2{i,5});
    DataMatrix(i,12)=mean(LengthDistribution2{i,6});
    
    DataMatrix(i,13)=std(LengthDistribution2{i,1});
    DataMatrix(i,14)=std(LengthDistribution2{i,2});
    DataMatrix(i,15)=std(LengthDistribution2{i,3});
    DataMatrix(i,16)=std(LengthDistribution2{i,4});
    DataMatrix(i,17)=std(LengthDistribution2{i,5});
    DataMatrix(i,18)=std(LengthDistribution2{i,6});
    
end

%% Plot circadian
% for i=1:nE
for i=([2,5,6,8]+20)
    
    figure
    subplot(2,1,1)
    for ii=1:nstateNames
        plot(circadianTime{i}/3600-lagVec(i),circadianStatistics{i}(:,ii),'Color',c(ii,:)),hold on
    end
    legend(stateNames)
    xtickVal=0:6:4*24;
    xticks(xtickVal)
    xticklabels(mod(xtickVal,24))
    xlim([circadianTime{i}(1),circadianTime{i}(end)]/3600-lagVec(i))
    ylim([0,100])
    grid on
    title(E(i).name(1:end-13),'interpreter','none')
    xlabel('Time (hours)')
    ylabel('Percentage')

    subplot(2,1,2)
    for ii=1:nstateNames
        plot(circadianFiltTime{i}/3600-lagVec(i),circadianFiltCell{i}(:,ii),'Color',c(ii,:)),hold on
    end
    legend(stateNames,'autoupdate','off')
    xtickVal=0:6:4*24;
    xticks(xtickVal)
    xticklabels(mod(xtickVal,24))
    xlim([circadianFiltTime{i}(1),circadianFiltTime{i}(end)]/3600-lagVec(i))
    grid on
    title([E(i).name(1:end-13),' -Lowpass filtered'],'interpreter','none')
    N=1e3;
    ylim([0,100])
    YLIM=get(gca,'YLim');
    contourT=linspace(-11,109,1e3);
    contourP=YLIM(2)*ones(1,N);
    contourZ=zeros(1,N);
    colormap(gray)
    contourC=repmat([linspace(1,0,N/10),linspace(0,1,N/10)],[1,5]);
    surface([contourT;contourT],[contourP;contourP],[contourZ;contourZ],[contourC;contourC],...
            'facecol','no',...
            'edgecol','interp',...
            'linew',5);
    xlabel('Time (hours)')
    ylabel('Percentage')

end
%%
figure
for iii=1:nstateNames
subplot(1,3,iii)
startPlot=zeros(nE,1);
endPlot=zeros(nE,1);
for i=1:nE
% for i=[1,2,3,4,7,8,24,25]
    set(gcf,'name',num2str(cut_off))
    plot(circadianFiltTime{i}/3600-lagVec(i),circadianFiltCell{i}(:,iii),'color',c(strcmp(DataMatrix_Genotype{i},genotypeNames),:)), hold on;
    startPlot(i)=circadianFiltTime{i}(1)/3600  -lagVec(i);
    endPlot(i)=  circadianFiltTime{i}(end)/3600-lagVec(i);
end
xtickVal=0:6:4*24;
xticks(xtickVal)
xticklabels(mod(xtickVal,24))
xlim([min(startPlot),max(endPlot)])
XLIM=get(gca,'XLim');
ylim([0,100])
YLIM=get(gca,'YLim');
N=1e3;
contourT=linspace(-11,109,1e3);
contourP=YLIM(2)*ones(1,N);
contourZ=zeros(1,N);
colormap(gray)
contourC=repmat([linspace(1,0,N/10),linspace(0,1,N/10)],[1,5]);
surface([contourT;contourT],[contourP;contourP],[contourZ;contourZ],[contourC;contourC],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',5);
    xlabel('Time (hours)')
    grid on
end
legendOb=gobjects(ngenotypeNames,1);
for iii=1:nstateNames
    subplot(1,3,iii), hold on
    for g=1:ngenotypeNames
        legendOb(g)=plot(NaN,'color',c(g,:));
    end
    legend(legendOb,genotypeNames,'autoupdate','off')
end
subplot(1,3,1)
ylabel('Awake (%)')
title('Awake as a percentage of segment time')
subplot(1,3,2)
ylabel('NREM (%)')
title('NREM as a percentage of segment time')
subplot(1,3,3)
ylabel('REM (%)')
title('REM as a percentage of sleep')

% Check locations of first local minimum and maximimum -Light=[7,19] is okay
locMin=zeros(nE,nstateNames);
locMax=zeros(nE,nstateNames);
circadianStartTime=zeros(nE,1);
circadianEndTime=zeros(nE,1);
for i=1:nE
    for iii=1:nstateNames
        if sum(islocalmin(circadianFiltCell{i}(:,iii)))~=0
        locMin(i,iii)=circadianFiltTime{i}(find(islocalmin(circadianFiltCell{i}(:,iii)),1,'first'))/3600-lagVec(i);
        else
        locMin(i,iii)=NaN;   
        end
        if sum(islocalmax(circadianFiltCell{i}(:,iii)))~=0
        locMax(i,iii)=circadianFiltTime{i}(find(islocalmax(circadianFiltCell{i}(:,iii)),1,'first'))/3600-lagVec(i);
        else
        locMax(i,iii)=NaN;   
        end
    end
    
    circadianStartTime(i)=circadianFiltTime{i}(1);
    circadianEndTime  (i)=circadianFiltTime{i}(end);
    
end
%%
figure
SEM_transparency=0.2;
circadianMeanSEM =cell(ngenotypeNames,1);
circadianMeanMean=cell(ngenotypeNames,1);
circadianMeanTime=cell(ngenotypeNames,1);
circadianFiltTime2=cell(nE,1);
for i=1:nE
    circadianFiltTime2{i}=circadianFiltTime{i}-lagVec(i)*3600;
end
for g=1:ngenotypeNames
    
    genLog=strcmp(DataMatrix_Genotype,genotypeNames{g});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     genLog=genLog&strcmp(DataMatrix_Day,'D1');
%     genLog=genLog&strcmp(DataMatrix_Day,'E1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    circadianMeanTime{g}=linspace(min(circadianStartTime(genLog)-lagVec(genLog)*3600),max(circadianEndTime(genLog)-lagVec(genLog)*3600),floor((max(circadianEndTime(genLog)-lagVec(genLog)*3600)-min(circadianStartTime(genLog)-lagVec(genLog)*3600))*fs_circadian))';
    nCircMean=length(circadianMeanTime{g});
    circadianMeanDT=circadianMeanTime{g}(2)-circadianMeanTime{g}(1);
    circadianFiltList=cat(1,circadianFiltCell{genLog});
    circadianTimeList=cat(1,circadianFiltTime2{genLog});
    circadianMeanIndex=round((circadianTimeList-circadianMeanTime{g}(1))/circadianMeanDT+1);

    circadianMeanNumel=accumarray(circadianMeanIndex,1,[nCircMean,1]);
    circadianMeanNumelSqrt=sqrt(circadianMeanNumel);
    circadianMeanSEM {g}=zeros(nCircMean,nstateNames);
    circadianMeanMean{g}=zeros(nCircMean,nstateNames);
    for iii=1:nstateNames
        circadianMeanSum  =accumarray(circadianMeanIndex,circadianFiltList(:,iii)',[nCircMean,1]);
        circadianMeanSTD  =accumarray(circadianMeanIndex,circadianFiltList(:,iii)',[nCircMean,1],@std);
         circadianMeanSEM{g}(:,iii)=circadianMeanSTD./circadianMeanNumelSqrt;
        circadianMeanMean{g}(:,iii)=circadianMeanSum./circadianMeanNumel;
        
        % filter Mean and SEM
         circadianMeanSEM{g}(:,iii)=LowPassfilter(2, cut_off, fs_circadian,  circadianMeanSEM{g}(:,iii));
        circadianMeanMean{g}(:,iii)=LowPassfilter(2, cut_off, fs_circadian, circadianMeanMean{g}(:,iii));
        
        subplot(1,nstateNames,iii),hold on
        plot(circadianMeanTime{g}/3600,circadianMeanMean{g}(:,iii),'color',c(g,:))
        patch('XData',[circadianMeanTime{g}/3600;...
                       circadianMeanTime{g}(end:-1:1)/3600;...
                       circadianMeanTime{g}(1)/3600],...
              'YData',[circadianMeanMean{g}(:,iii)-circadianMeanSEM{g}(:,iii);...
                       circadianMeanMean{g}(end:-1:1,iii)+circadianMeanSEM{g}(end:-1:1,iii);...
                       circadianMeanMean{g}(1,iii)-circadianMeanSEM{g}(1,iii)],...
           'FaceColor',c(g,:),...
           'EdgeColor',c(g,:),...
           'FaceVertexCData',c(g,:),...
           'FaceAlpha',SEM_transparency,...
           'EdgeAlpha',SEM_transparency,...
           'FaceVertexAlphaData',SEM_transparency)
    end
    
end
legendOb=gobjects(ngenotypeNames,1);
for iii=1:nstateNames
    subplot(1,nstateNames,iii)
    for g=1:ngenotypeNames
        legendOb(g)=plot(NaN,NaN,'Color',c(g,:));
    end
    legend(legendOb,genotypeNames,'autoupdate','off')
    
    xtickVal=0:6:4*24;
    xticks(xtickVal)
    xticklabels(mod(xtickVal,24))
    xlim([min(circadianStartTime),max(circadianEndTime)]/3600)
    XLIM=get(gca,'XLim');
%     ylim([0,100])
    ylim('auto')
    YLIM=get(gca,'YLim');
    N=1e3;
    contourT=linspace(-11,109,1e3);
    contourP=YLIM(2)*ones(1,N);
    contourZ=zeros(1,N);
    colormap(gray)
    contourC=repmat([linspace(1,0,N/10),linspace(0,1,N/10)],[1,5]);
    surface([contourT;contourT],[contourP;contourP],[contourZ;contourZ],[contourC;contourC],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',5);
    xlabel('Time (hours)')
    grid on
end
subplot(1,3,1)
ylabel('Awake (%)')
title('Awake as a percentage of segment time')
subplot(1,3,2)
ylabel('NREM (%)')
title('NREM as a percentage of segment time')
subplot(1,3,3)
ylabel('REM (%)')
title('REM as a percentage of sleep')

%% Plot boxplot

% Plot all 18 histograms
H=cell(18,1);
P=cell(18,1);
P2=cell(18,1);
P_Significance=zeros(18,1);
P_Significance_anova=zeros(18,1);
CI=cell(18,1);
STATS=cell(18,1);
p_threshold=0.05;

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
      if nowtrescue
          Gen={'wt';'hetrescue';'het'};
          nGen=length(Gen);
boxplot([DataMatrix(strcmp(DataMatrix_Genotype2,Gen{1}),i);...
              DataMatrix(strcmp(DataMatrix_Genotype2,Gen{2}),i);...
              DataMatrix(strcmp(DataMatrix_Genotype2,Gen{3}),i)],...
              [repmat({['WT (n=',num2str(sum(strcmp(DataMatrix_Genotype2,Gen{1}))),')']},[sum(strcmp(DataMatrix_Genotype2,Gen{1})),1]);...
               repmat({['het+ (n=',num2str(sum(strcmp(DataMatrix_Genotype2,Gen{2}))),')']},[sum(strcmp(DataMatrix_Genotype2,Gen{2})),1]);...
               repmat({['het- (n=',num2str(sum(strcmp(DataMatrix_Genotype2,Gen{3}))),')']},[sum(strcmp(DataMatrix_Genotype2,Gen{3})),1])],...
               'Colors',c([2,3,1],:))
      else
          Gen={'wt';'wtrescue';'hetrescue';'het'};
          nGen=length(Gen);
      boxplot([DataMatrix(strcmp(DataMatrix_Genotype2,Gen{1}),i);...
               DataMatrix(strcmp(DataMatrix_Genotype2,Gen{2}),i);...
               DataMatrix(strcmp(DataMatrix_Genotype2,Gen{3}),i);...
               DataMatrix(strcmp(DataMatrix_Genotype2,Gen{4}),i)],...
              [repmat({['WT (n=',num2str(sum(strcmp(DataMatrix_Genotype2,Gen{1}))),')']},[sum(strcmp(DataMatrix_Genotype2,Gen{1})),1]);...
               repmat({['WT+ (n=',num2str(sum(strcmp(DataMatrix_Genotype2,Gen{2}))),')']},[sum(strcmp(DataMatrix_Genotype2,Gen{2})),1]);...
               repmat({['het+ (n=',num2str(sum(strcmp(DataMatrix_Genotype2,Gen{3}))),')']},[sum(strcmp(DataMatrix_Genotype2,Gen{3})),1]);...
               repmat({['het- (n=',num2str(sum(strcmp(DataMatrix_Genotype2,Gen{4}))),')']},[sum(strcmp(DataMatrix_Genotype2,Gen{4})),1])],...
               'Colors',c([2,4,3,1],:))
      end
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
    
    if nowtrescue
    Data_Current={DataMatrix(strcmp(DataMatrix_Genotype2,Gen{1}),i);...
                  DataMatrix(strcmp(DataMatrix_Genotype2,Gen{2}),i);...
                  DataMatrix(strcmp(DataMatrix_Genotype2,Gen{3}),i)};
    else
    Data_Current={DataMatrix(strcmp(DataMatrix_Genotype2,Gen{1}),i);...
                  DataMatrix(strcmp(DataMatrix_Genotype2,Gen{2}),i);...
                  DataMatrix(strcmp(DataMatrix_Genotype2,Gen{3}),i);...
                  DataMatrix(strcmp(DataMatrix_Genotype2,Gen{4}),i)};
    end
              
    for g=1:length(Data_Current)
    for gg=1:length(Data_Current)
        [H{i}(g,gg),P{i}(g,gg),CI{i}{g,gg},STATS{i}{g,gg}]=ttest2(Data_Current{g},Data_Current{gg},'vartype','unequal');
    end
    end
    P_Significance(i)=round(sum(sum(P{i}<p_threshold))/2);
    
% anova start
anova_data=Data_Current;
g1=cell(nGen,1);
g2=cell(nGen,1);

for g=1:nGen
    if strcmp(Gen{g},'wt')||strcmp(Gen{g},'wtrescue')
        g1{g}=repmat({'wt'},[1,length(anova_data{g})]);
    else
        g1{g}=repmat({'het'},[1,length(anova_data{g})]);
    end
    if strcmp(Gen{g},'wtrescue')||strcmp(Gen{g},'hetrescue')
        g2{g}=repmat({'yes'},[1,length(anova_data{g})]);
    else
        g2{g}=repmat({'no'},[1,length(anova_data{g})]);
    end
end
anova_data=cat(1,anova_data{:})';
g1=[g1{:}];
g2=[g2{:}];

[~,~,stats,~]=anovan(anova_data,{g1,g2},'model','interaction','varnames',{'genotype','rescue'},'display','off');
[results,~,~,gnames] = multcompare(stats,'Dimension',[1 2],'display','off'); 

for g=1:size(results,2)
    
    if strcmp(gnames(results(g,1)),{'genotype=wt,rescue=no' })
        gname_temp1='wt';
    elseif strcmp(gnames(results(g,1)),{'genotype=wt,rescue=yes' })
        gname_temp1='wtrescue';
    elseif strcmp(gnames(results(g,1)),{'genotype=het,rescue=yes' })
        gname_temp1='hetrescue';
    elseif strcmp(gnames(results(g,1)),{'genotype=het,rescue=no' })
        gname_temp1='het';
    end
    
    if strcmp(gnames(results(g,2)),{'genotype=wt,rescue=no' })
        gname_temp2='wt';
    elseif strcmp(gnames(results(g,2)),{'genotype=wt,rescue=yes' })
        gname_temp2='wtrescue';
    elseif strcmp(gnames(results(g,2)),{'genotype=het,rescue=yes' })
        gname_temp2='hetrescue';
    elseif strcmp(gnames(results(g,2)),{'genotype=het,rescue=no' })
        gname_temp2='het';
    end
    
    P2{i}(find(strcmp(Gen,gname_temp1)),find(strcmp(Gen,gname_temp2)))=results(g,end);
    P2{i}(find(strcmp(Gen,gname_temp2)),find(strcmp(Gen,gname_temp1)))=results(g,end);
    
end
P2{i}(eye(nGen)==1)=1;
% anova end
P_Significance_anova(i)=round(sum(sum(P2{i}<p_threshold))/2);
    
    
end
P_Significance=[P_Significance,(1:18)'];
P_Significance_anova=[P_Significance_anova,(1:18)'];

%% Save in table
% Average circadian rhythms with grid (and use distance)


% DataMatrix(:,4:6)=[(DataMatrix(:,4:5)./repmat(sum(DataMatrix(:,4:6),2),[1,2]))*100,(DataMatrix(:,6)./sum(DataMatrix(i,5:6)))*100];
% DataMatrix(:,1:3)=[(DataMatrix(:,1:2)./repmat(sum(DataMatrix(:,1:3),2),[1,2]))*100,(DataMatrix(:,3)./sum(DataMatrix(i,2:3)))*100];