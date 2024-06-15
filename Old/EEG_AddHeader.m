%% Add header
MouseDir=dir('C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Analyzed\*_*_*_Analyzed.mat');
nMouse=length(MouseDir);
check=false(nMouse,1);

for i=1:nMouse

load(['C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Clean\',...
MouseDir(i).name(1:end-13),'.mat'], 'header')

HEADER1=header;

save(['C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Analyzed\',...
MouseDir(i).name(1:end-13),'_Analyzed.mat'],'header','-append')

load(['C:\Users\enzo\Downloads\Study\Current Courses\BEP\Data\EEG\Data_Analyzed\',...
MouseDir(i).name(1:end-13),'_Analyzed.mat'], 'header')

HEADER2=header;

check(i)=strcmp(HEADER1.starttime,HEADER2.starttime);

end
