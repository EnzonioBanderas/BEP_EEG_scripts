% Marcel van Velze 19.09.2017
% marcelbvv@gmail.com
% Enzo Nio
% enzonio@hotmail.com

% .edf converter
%
% Info: This script can be used to convert .edf files converted from 
%       .hex files by Neurologger-II software to .mat files.
%       Place .edf files in one folder and select that folder.
clearvars, close all
%% Go to and get experiment folder directory
experiment_folder = uigetdir('','Select folder containing .edf file');
cd (experiment_folder)
experiments = dir('*.edf');

%% Convert .edf to .mat 
bar = waitbar(0, 'Converting Data, please wait');
for i = 1:length(experiments)

   waitbar((i-1)/length(experiments), bar, 'Converting Data, please wait');
   
   % get name of current experiment in loop
   experiment_name=experiments(i).name;
   
   % convert with edfread
   [header, data] = edfread(experiment_name);
   
   % transpose of matrix (MATLAB uses column major order)
   data=data';
   
   % Seperate data into a time x type matrix for data with different
   % sampling rates
   Electrical.CH1234=data(:,1:4);
   REF12.data=data(:,5:6);
   IR.data=data(:,7);
   Acceleration.XYZ=data(:,8:10);
   Unknown.data=data(:,11);
   
   % Add sampling rates to data structure
   Electrical.fs=header.frequency(1);
   REF12.fs=header.frequency(5);
   IR.fs=header.frequency(7);
   Acceleration.fs=header.frequency(8);
   Unknown.fs=header.frequency(11);
   
   % truncation of data with a lower sampling rate than electric data 
   % (other data than electric is padded with zeros at the end)
   T=size(Electrical.CH1234,1)/Electrical.fs;
   trunc_IR=round(T*IR.fs);
   trunc_Acceleration=round(T*Acceleration.fs);
   trunc_Unknown=round(T*Unknown.fs);
   
   IR.data_bla=IR.data(1:trunc_IR);
   Acceleration.XYZ=Acceleration.XYZ(1:trunc_Acceleration,:);
   Unknown.data=IR.data(1:trunc_Unknown);
   
   % Save data to .mat file with the same name
   save(experiment_name(1:end-4),'Electrical','REF12','IR','Acceleration','Unknown','header','-v7.3')
   
   % Clear all unnecessary variables to save memory
   clearvars -except time_convert i experiments bar
   
end
close(bar)