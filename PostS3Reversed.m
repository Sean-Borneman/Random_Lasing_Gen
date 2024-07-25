% HOUSEKEEPING
close all;
clearvars ;

fontSize = 20;
%% GET THE DATA
dataPath = 'C:/Users/Sean/Documents/Programing/Python/PURDUE/ControlledRunning/Power_238[uJcm2]/';
flatLasingQ = table2array(readtable([dataPath 'flatLasngQ[uJcm2].csv'],'ReadVariableNames', false)); %this data includes thye diagnal


matrixB = abs(flatLasingQ);

figure('Position',[50 165 800 550]);
histogram(matrixB,'EdgeColor','k','EdgeAlpha',0.3,'FaceColor','m');
title('$\textrm{Histogram } q_{lasing}$','Interpreter','latex','FontSize',fontSize);
xlabel('$q$','Interpreter','latex','FontSize',fontSize)
ylabel('$Counts$','Interpreter','latex','FontSize',fontSize)


lasngPd = fitdist(flatLasingQ,'Kernel','Kernel','normal'); %Fit the flat lasing data to a distribution

%% Unflatten
a =100;% numberOfReplica
lasingQ = reshape(flatLasingQ, a,a);

%% Step 1: Extract the Delta
lasngDelta = sqrtm(lasingQ); %the noraml delta is 100x158 not 100x100
                              %effectively this cahnges the lasing wavelengths 

%% Step 2: 
%now you have to add the mean which seems like a lot of info to suddenly
%add. maybe add the real mean??
