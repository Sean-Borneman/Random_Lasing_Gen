%% Clean UP
close all;
clear all;
clc;
%% Load Data
load('SpectralStatData.mat');
POWER = 8;
Data.power(POWER)*1000

%% PLOT single pulse
%plot(Data.Wavel, Data.singleSpec(1).data(:,1))
%plot(Data.Wavel ,Data.singleSpec(7).data)
imagesc(Data.singleSpec(POWER).data)
%writematrix(Data.Wavel,"Frequencies.csv", 'Delimiter',',')
%FILENAME = "RawData"+ round(Data.power(POWER)*1000)+ "[uJcm2].csv";%
%writematrix(Data.singleSpec(POWER).data,FILENAME,'Delimiter',',');