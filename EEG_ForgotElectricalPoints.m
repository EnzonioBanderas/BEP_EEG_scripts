clearvars, close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EF = uigetdir('','Select Experiment Folder (EF)');
cd(EF)
E=dir('*_Analyzed.mat');
nE=length(E);

EleAccMultiple=5;

for i=1:nE
%% Take maximum of probabilities to create state logicals
load(E(i).name,'Acceleration','Electrical','Electrical_P')
Electrical.points_FA=logical2points(~isnan(Electrical_P(:,1)));
[~,P_FA_index]=max(Electrical_P,[],2);
P_FA_index(isnan(Electrical_P(:,1)))=NaN;
Electrical.points_AA=logical2points(P_FA_index==1);
Electrical.points_NREM=logical2points(P_FA_index==2);
Electrical.points_REM=logical2points(P_FA_index==3);
Electrical.points_QA=logical2points(P_FA_index==4);
clearvars P_FA_index

% Adjust Electrical points (probably has no effect for standard Settings)
Electrical.points_FA(:,1)     =Electrical.points_FA(:,1)-1; Electrical.points_FA=round(Electrical.points_FA/EleAccMultiple)*EleAccMultiple; Electrical.points_FA(:,1)=Electrical.points_FA(:,1)+1;
Electrical.points_AA(:,1)     =Electrical.points_AA(:,1)-1; Electrical.points_AA=round(Electrical.points_AA/EleAccMultiple)*EleAccMultiple; Electrical.points_AA(:,1)=Electrical.points_AA(:,1)+1;
Electrical.points_NREM(:,1)   =Electrical.points_NREM(:,1)-1; Electrical.points_NREM=round(Electrical.points_NREM/EleAccMultiple)*EleAccMultiple; Electrical.points_NREM(:,1)=Electrical.points_NREM(:,1)+1;
Electrical.points_REM(:,1)    =Electrical.points_REM(:,1)-1; Electrical.points_REM=round(Electrical.points_REM/EleAccMultiple)*EleAccMultiple; Electrical.points_REM(:,1)=Electrical.points_REM(:,1)+1;
Electrical.points_QA(:,1)     =Electrical.points_QA(:,1)-1; Electrical.points_QA=round(Electrical.points_QA/EleAccMultiple)*EleAccMultiple; Electrical.points_QA(:,1)=Electrical.points_QA(:,1)+1;

% Convert Electrical points to Acceleration points
Acceleration.points_FA     =[Electrical.points_FA(:,1)+EleAccMultiple-1,Electrical.points_FA(:,2)]/EleAccMultiple;
Acceleration.points_AA     =[Electrical.points_AA(:,1)+EleAccMultiple-1,Electrical.points_AA(:,2)]/EleAccMultiple;
Acceleration.points_NREM   =[Electrical.points_NREM(:,1)+EleAccMultiple-1,Electrical.points_NREM(:,2)]/EleAccMultiple;
Acceleration.points_REM    =[Electrical.points_REM(:,1)+EleAccMultiple-1,Electrical.points_REM(:,2)]/EleAccMultiple;
Acceleration.points_QA     =[Electrical.points_QA(:,1)+EleAccMultiple-1,Electrical.points_QA(:,2)]/EleAccMultiple;

clearvars 'Electrical_P'

save(E(i).name,'Electrical','Acceleration','-append')   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end