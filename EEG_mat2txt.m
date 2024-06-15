clearvars, close all
% Get files and define settings
cd('C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data')
[FileName,PathName]=uigetfile('*_*_*.mat','Select .mat Neurologger file');
cd(PathName)
load(FileName,'Acceleration','Electrical')
Mouse=FileName(1:end-4);
Settings.HP_acc = 1; % Hz
Settings.order_acc = 2; % order of filter
Settings.nSplit=4; % split in nSplit parts
Settings.MaxTime=12; % hours
EleAccMultiple=round(Electrical.fs/Acceleration.fs);

% Calculate Dynamic Acceleration vector length
Acceleration.dyn=[HighPassfilter(Settings.order_acc,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,1)),...
HighPassfilter(Settings.order_acc,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,2)),...
HighPassfilter(Settings.order_acc,Settings.HP_acc,Acceleration.fs,Acceleration.XYZ(:,3))];
Acceleration.dyn=vecnorm(Acceleration.dyn,2,2);

% Rejection of data for 20s periods of EEG that come within 10 \muV of maximum
% selected=~conv(abs(Electrical.CH1234(:,2))>max(Electrical.CH1234(:,2))-10,true(round(20*Electrical.fs),1),'same');
Right_M=Electrical.CH1234(:,2); % choose right-M
DynAccLength=Acceleration.dyn;

% Number of points for split segments
nSegment=size(Electrical.points,1);
nPoint=floor((round(Settings.MaxTime*60*60*Electrical.fs))/5)*5;
nSplit=floor((Electrical.points(:,2)-Electrical.points(:,1)+1)/nPoint);

count=0;
for i=1:nSegment
    
for ii=1:nSplit(i)+1
count=count+1;
    
    
% Define or update index
if ii==1
    
index_start=Electrical.points(i,1); 
index_end=Electrical.points(i,1)+nPoint-1;

else
index_start=index_start+nPoint; index_end=index_end+nPoint;

end

if index_end>Electrical.points(i,2)
    index_end=Electrical.points(i,2);
end

index=index_start:index_end;
EEG=Right_M(index); % choose right-M


index2=(index_start+EleAccMultiple-1)/5:index_end/5;
EMG=resample(DynAccLength(index2),EleAccMultiple,1); % fake EMG

% Write EEG and 'EMG' to .txt files
cd('C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Kreuzer')
dlmwrite([Mouse,'_',num2str(index_start),'-',num2str(index_end),'_EEG',num2str(count),'.txt'],EEG)
dlmwrite([Mouse,'_',num2str(index_start),'-',num2str(index_end),'_EMG',num2str(count),'.txt'],EMG)

end

end