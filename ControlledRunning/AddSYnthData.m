close all
clear all
clc
%% GET THE DATA
dataPath = 'C:/Users/Sean/Documents/Programing/Python/PURDUE/ControlledRunning/';
%dataPath = 'C:/Users/Sean/Documents/Programing/Python/PURDUE/RealData/';

load([dataPath 'SpectralStatData.mat']);
%% add 
freqMix = table2array(readtable('freqMix.csv'));
Data.singleSpec(1).data = freqMix
plot(Data.Wavel ,Data.singleSpec(1).data)
