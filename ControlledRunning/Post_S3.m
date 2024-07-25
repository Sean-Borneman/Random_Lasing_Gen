%% POST-PROCESSING OF THE RANDOM LASING DATA
% DATA STRUCTURE
% Data.LaserEn, 
% Data.LaserD, 
% Data.LaserT, 
% Data.Reprate, 
% Data.NDfilter, 
% Data.VarND, 
% Data.Slitwidth, 
% Data.Exposuretime, 
% Data.Grating, 
% Data.Wavel, 
% Data.Path, 
% Data.m, 
% Data.th0, 
% Data.spec, 
% Data.power, 
% Data.angle, 
% Data.singleSpec, 
% Data.ADCrate, 
% Data.totSpec

% HOUSEKEEPING
close all;
%clearvars ;

%% GET THE DATA
dataPath = 'C:/Users/Sean/Documents/Programing/Python/PURDUE/ControlledRunning/';
%dataPath = 'C:/Users/Sean/Documents/Programing/Python/PURDUE/RealData/';

load([dataPath 'SpectralStatData.mat']);
% 
freqMix =table2array(readtable('freqMIX238.csv'));% table2array(readtable('avrgwDist.csv'));
%Data.singleSpec(1).data = freqMix
%% CONSTANTS
close all 
binsNumber = 30;
fontSize   = 20;
wavelength = Data.Wavel;
lasngIndex = (find(wavelength<=550,1,'last'):find(wavelength>=557,1)); % 7 nm centered around lasing wavel
spontIndex = (find(wavelength<=536,1,'last'):find(wavelength>=543,1)); % 7 nm centered around 
powerIndex = (1:length(Data.power));
lasngCentr = [556,555,554.5,554,553.5,553,552.5,552,551.5,551,551,551,551,551];

%% MAIN LOOP - SCAN ALL THE 

numberOfReplica = size(Data.singleSpec(1).data,2);
numberOfWavelen = size(Data.singleSpec(1).data,1);
numberOfDataPow = length(Data.power); 

sizeQ = numberOfDataPow*numberOfDataPow-numberOfDataPow;

flatSpontQ     = zeros(sizeQ,1);
flatLasngQ     = zeros(sizeQ,1);
flatNormSpontQ = zeros(sizeQ,1);
flatNormLasngQ = zeros(sizeQ,1);


for p = 1:1%numberOfDataPow

close all;

numPOWER = round(Data.power(p)*1000);
strPOWER = num2str(numPOWER);

dataFolder = ['Power_' strPOWER '[uJcm2]']; 

mkdir(dataFolder);

%% STEP 1: Take the 7-nm band centered around the lasing wavelength
lasngIndex = (find(wavelength<=(lasngCentr(p)-3.5),1,'last'):find(wavelength>=(lasngCentr(p)+3.5),1));

%% STEP 2: GET AND STORE SPECTROGAMM MATRIX Ip FOR A GIVEN FLUENCE P  
Ip      = Data.singleSpec(p).data;
spontIp = Ip(spontIndex,:);
lasngIp = Ip(lasngIndex,:);

FILENAME = [dataPath dataFolder '/spontIp_' num2str(Data.power(p)) '[uJcm2].csv'];
writematrix(spontIp,FILENAME,'Delimiter',',');
FILENAME = [dataPath dataFolder '/lasngIp_' num2str(Data.power(p)) '[uJcm2].csv'];
writematrix(lasngIp,FILENAME,'Delimiter',',');

%% STEP 3: AVERAGE ALL REPLICAS IN Ip TO GET A MEAN VECTOR, 
%% REPLICATED INTO A MATRIX Ibar AND STORE IT IN A FILE
spontIbar = repmat(mean(spontIp,2),[1,size(spontIp,2)]);
lasngIbar = repmat(mean(lasngIp,2),[1,size(lasngIp,2)]);

FILENAME  = [dataPath dataFolder '/spontIp_' num2str(Data.power(p)) '[uJcm2].csv'];
writematrix(spontIp,FILENAME,'Delimiter',',');
FILENAME  = [dataPath dataFolder '/lasngIp_' num2str(Data.power(p)) '[uJcm2].csv'];
writematrix(lasngIp,FILENAME,'Delimiter',',');

%% STEP 4: GET THE DIFFERENCE MATRICES \Delta FOR A GIVEN FLUENCE P 
spontDelta     = spontIp - spontIbar;
lasngDelta     = lasngIp - lasngIbar;

FILENAME = [dataPath dataFolder '/spontDelta_' num2str(Data.power(p)) '[uJcm2].csv'];
writematrix(spontDelta,FILENAME,'Delimiter',',');
FILENAME = [dataPath dataFolder '/lasngDelta_' num2str(Data.power(p)) '[uJcm2].csv'];
writematrix(lasngDelta,FILENAME,'Delimiter',',');

%% STEP 5: NORMALIZE THE DIFFERENCE MATRICES \Delta 
wavelengthNumber     = size(spontDelta,1);
absSpontDelta        = repmat(sqrt(sum(spontDelta.^2,1)),[wavelengthNumber,1]);
normalizedSpontDelta = spontDelta./absSpontDelta; 

wavelengthNumber     = size(lasngDelta,1);
absLasngDelta        = repmat(sqrt(sum(lasngDelta.^2,1)),[wavelengthNumber,1]);
normalizedLasngDelta = lasngDelta./absLasngDelta;

spontQ = abs(spontDelta'*spontDelta);
lasngQ = abs(lasngDelta'*lasngDelta);

normalizedSpontQ = abs(normalizedSpontDelta'*normalizedSpontDelta);
normalizedLasngQ = abs(normalizedLasngDelta'*normalizedLasngDelta);


%% IGNORE THE DATA AT THE DIAGONAL
%  AND FLATTEN THE ARRAY


a  = numberOfReplica;
a1 = a+1;
a2 = a*a;


flatSpontQ = reshape(spontQ, a2,1);
flatLasngQ = reshape(lasngQ, a2,1);

flatNormSpontQ = reshape(normalizedSpontQ,a2,1);
flatNormLasngQ = reshape(normalizedLasngQ,a2,1);
          

% flatSpontQ(1:a1:a2)=[];
% flatLasngQ(1:a1:a2)=[];
% 
% flatNormSpontQ(1:a1:a2)=[];
% flatNormLasngQ(1:a1:a2)=[];


FILENAME = [dataPath dataFolder '/flatSpontQ_' num2str(Data.power(p)) '[uJcm2].csv'];
writematrix(flatSpontQ,FILENAME,'Delimiter',',');
FILENAME = [dataPath dataFolder '/flatLasngQ_' num2str(Data.power(p)) '[uJcm2].csv'];
writematrix(flatLasngQ,FILENAME,'Delimiter',',');

FILENAME = [dataPath dataFolder '/flatNormSpontQ_' num2str(Data.power(p)) '[uJcm2].csv'];
writematrix(flatNormSpontQ,FILENAME,'Delimiter',',');
FILENAME = [dataPath dataFolder '/flatNormLasngQ_' num2str(Data.power(p)) '[uJcm2].csv'];
writematrix(flatNormLasngQ,FILENAME,'Delimiter',',');



%% first quadrant
figure;

subplot(1,5,(1:4));
imagesc(1:size(Ip,2),wavelength,Ip)
hold on
% x = [0, 101, 101, 0, 0];
% y = [wavelength(lasngIndex(1)),wavelength(lasngIndex(1)),wavelength(lasngIndex(end)),wavelength(lasngIndex(end)),wavelength(lasngIndex(1))];
% fill(x,y,'m','FaceAlpha',0.3)
yline(wavelength(lasngIndex(1)),'m','Linewidth',3)
yline(wavelength(lasngIndex(end)),'m','Linewidth',3)
yline(wavelength(spontIndex(1)),'r','Linewidth',3)
yline(wavelength(spontIndex(end)),'r','Linewidth',3)

title('$I^{(\alpha)}(\lambda_k)$','Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\lambda \textrm{[nm]}$','Interpreter','latex')
set(gca,'FontSize',fontSize)
pbaspect(gca,[1,1,1]);

Ibar = repmat(mean(Ip,2),[1,2]);

subplot(1,5,5)
hold on
imagesc(1:size(Ip,2),wavelength,Ibar);
pbaspect(gca,[0.1,2,0.1]);
title('$\bar{I}(\lambda)$','Interpreter','latex');
colorbar;

x = [0, 101, 101, 0, 0];
y = [wavelength(lasngIndex(1)),wavelength(lasngIndex(1)),wavelength(lasngIndex(end)),wavelength(lasngIndex(end)),wavelength(lasngIndex(1))];

yline(wavelength(lasngIndex(1)),'m','Linewidth',3)
yline(wavelength(lasngIndex(end)),'m','Linewidth',3)
yline(wavelength(spontIndex(1)),'r','Linewidth',3)
yline(wavelength(spontIndex(end)),'r','Linewidth',3)

set(gca,'YDir','reverse');
set(gca,'YLim',[min(wavelength),max(wavelength)]);
%set(gca,'XLim',[0.5*min(Ibar(:,1)),1.05*max(Ibar(:,1))]);
set(gca,'FontSize',fontSize)
set(gca,'xTick',[]);
set(gca,'yTick',[]);
ylabel('$\lambda \textrm{[nm]}$','Interpreter','latex');

FIG = [dataPath dataFolder '/I_pow' strPOWER]; 

savefig(gcf,[FIG '.fig']);
print(gcf,[FIG '.eps'],'-deps','-painters')
print(gcf,[FIG '.png'],'-dpng','-painters')


%% RAW DELTAS (w/o normalization)

matrixA = abs(spontDelta);
matrixB = abs(lasngDelta);

figure('Position',[50 165 550 785]);
subplot(2,1,1)
imagesc(1:100,wavelength(spontIndex),matrixA)
graymap = colormap('gray');
%magentaMap = graymap.*repmat([1,0,1],[256,1]); %commented because of strage matrix dimensions error
%redMap = graymap.*repmat([1,0,0],[256,1]);
%colormap(gca,redMap)
%colorbar;
title('$\Delta_{spont}^{(\alpha)}(\lambda_k)$','Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\lambda [nm]$','Interpreter','latex')
set(gca,'FontSize',fontSize)
pbaspect(gca,[1,1,1])

subplot(2,1,2)
imagesc(1:100,wavelength(lasngIndex),matrixB)
%colormap(gca,magentaMap) %comment because of abovce error
%colorbar;
title('$\Delta_{lasing}^{(\alpha)}(\lambda_k)$','Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\lambda [nm]$','Interpreter','latex')
set(gca,'FontSize',fontSize)
pbaspect(gca,[1,1,1])

FIG = [dataPath dataFolder '/rawDelta_pow' strPOWER]; 

savefig(gcf,[FIG '.fig']);
print(gcf,[FIG '.eps'],'-deps','-painters')
%print(gcf,[FIG '.png'],'-dpng','-painters')


%% NORMALIZED DELTAS

matrixA = abs(normalizedSpontDelta);
matrixB = abs(normalizedLasngDelta);

figure('Position',[50 165 550 785]);
subplot(2,1,1)
imagesc(1:100,wavelength(spontIndex),matrixA)
graymap = colormap('gray');
%magentaMap = graymap.*repmat([1,0,1],[256,1]);
%redMap = graymap.*repmat([1,0,0],[256,1]);
%colormap(gca,redMap)
%colorbar;
title('$\Delta_{spont}^{(\alpha)}(\lambda_k)$','Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\lambda [nm]$','Interpreter','latex')
set(gca,'FontSize',fontSize)
pbaspect(gca,[1,1,1])

subplot(2,1,2)
imagesc(1:100,wavelength(lasngIndex),matrixB)
%colormap(gca,magentaMap)
%colorbar;
title('$\Delta_{lasing}^{(\alpha)}(\lambda_k)$','Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\lambda [nm]$','Interpreter','latex')
set(gca,'FontSize',fontSize)
pbaspect(gca,[1,1,1])

FIG = [dataPath dataFolder '/normDelta_pow' strPOWER]; 

savefig(gcf,[FIG '.fig']);
print(gcf,[FIG '.eps'],'-deps','-painters')
print(gcf,[FIG '.png'],'-dpng','-painters')


%% RAW Q (w/o normalization)

matrixA = abs(spontQ);
matrixB = abs(lasngQ);

figure('Position',[50 165 550 785]);
subplot(2,1,1)
imagesc(matrixA);
title('$q_{spont}^{(\alpha\beta)}$','Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\beta$','Interpreter','latex')
%colormap(gca,redMap)
%colorbar;
pbaspect(gca,[1,1,1])
set(gca,'FontSize',fontSize)

subplot(2,1,2)
imagesc(matrixB);
title('$q_{lasing}^{(\alpha\beta)}$','Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\beta$','Interpreter','latex')
%colormap(gca,magentaMap)
%colorbar;
pbaspect(gca,[1,1,1])
set(gca,'FontSize',fontSize)

FIG = [dataPath dataFolder '/rawQ_pow' strPOWER]; 

savefig(gcf,[FIG '.fig']);
print(gcf,[FIG '.eps'],'-deps','-painters')
print(gcf,[FIG '.png'],'-dpng','-painters')


%% NORMALIZED Q

matrixA = abs(normalizedSpontQ);
matrixB = abs(normalizedLasngQ);

figure('Position',[50 165 550 785]);
subplot(2,1,1)
imagesc(matrixA);
title('$q_{spont}^{(\alpha\beta)}$','Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\beta$','Interpreter','latex')
%colormap(gca,redMap)
%colorbar;
pbaspect(gca,[1,1,1])
set(gca,'FontSize',fontSize)

subplot(2,1,2)
imagesc(matrixB);
title('$q_{lasing}^{(\alpha\beta)}$','Interpreter','latex');
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\beta$','Interpreter','latex')
%colormap(gca,magentaMap)
%colorbar;
pbaspect(gca,[1,1,1])
set(gca,'FontSize',fontSize)

FIG = [dataPath dataFolder '/normQ_pow' strPOWER]; 

savefig(gcf,[FIG '.fig']);
print(gcf,[FIG '.eps'],'-deps','-painters')
print(gcf,[FIG '.png'],'-dpng','-painters')


%% Raw Histograms
close all

matrixA = abs(flatSpontQ);
matrixB = abs(flatLasngQ);

%dbin = 1/binsNumber;
%edges = -dbin/2:dbin:(1+dbin/2);
figure('Position',[50 165 800 550]);
subplot(1,2,1)
histogram(matrixA,'EdgeColor','k','EdgeAlpha',0.3,'FaceColor','r');
title('$\textrm{Histogram } q_{spont}$','Interpreter','latex');
xlabel('$q$','Interpreter','latex')
ylabel('$Counts$','Interpreter','latex')
pbaspect(gca,[1,1,1])
set(gca,'FontSize',fontSize)

subplot(1,2,2)
histogram(matrixB,'EdgeColor','k','EdgeAlpha',0.3,'FaceColor','m');
lasingHistogram = matrixB;
title('$\textrm{Histogram } q_{lasing}$','Interpreter','latex','FontSize',fontSize);
xlabel('$q$','Interpreter','latex','FontSize',fontSize)
ylabel('$Counts$','Interpreter','latex','FontSize',fontSize)
%set(gca,'XLim',[-dbin/2,1+dbin/2])
pbaspect(gca,[1,1,1])
set(gca,'FontSize',fontSize)


FIG = [dataPath dataFolder '/rawHist_pow' strPOWER]; 

savefig(gcf,[FIG '.fig']);
print(gcf,[FIG '.eps'],'-deps','-painters')
print(gcf,[FIG '.png'],'-dpng','-painters')

%% NORMALIZED HISTOGRAM

matrixA = abs(flatNormSpontQ);
matrixB = abs(flatNormLasngQ);

%dbin = 1/binsNumber;
%edges = -dbin/2:dbin:(1+dbin/2);
figure('Position',[50 165 800 550]);
subplot(1,2,1)
histogram(matrixA,'EdgeColor','k','EdgeAlpha',0.3,'FaceColor','r');
title('$\textrm{Histogram } q_{spont}$','Interpreter','latex');
xlabel('$q$','Interpreter','latex')
ylabel('$Counts$','Interpreter','latex')
pbaspect(gca,[1,1,1])
set(gca,'FontSize',fontSize)
%set(gca,'XLim',[-dbin/2,1+dbin/2])

subplot(1,2,2)
histogram(matrixB,'EdgeColor','k','EdgeAlpha',0.3,'FaceColor','m');
title('$\textrm{Histogram } q_{lasing}$','Interpreter','latex','FontSize',fontSize);
xlabel('$q$','Interpreter','latex','FontSize',fontSize)
ylabel('$Counts$','Interpreter','latex','FontSize',fontSize)
%set(gca,'XLim',[-dbin/2,1+dbin/2])
pbaspect(gca,[1,1,1])
set(gca,'FontSize',fontSize)


FIG = [dataPath dataFolder '/normHist_pow' strPOWER]; 

savefig(gcf,[FIG '.fig']);
print(gcf,[FIG '.eps'],'-deps','-painters')
print(gcf,[FIG '.png'],'-dpng','-painters')



%% other Histogram
DATA = [flatLasngQ; -flatLasngQ];


histFigure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
% AXES
axesHistogram1 = axes('Parent',histFigure1); hold(axesHistogram1,'on');
histHeader1 = histogram(DATA,'EdgeColor','k','EdgeAlpha',0.3,'FaceColor','m');
title(['Laser Fluence ' strPOWER ' \muJ/cm^2'],'Interpreter','tex','FontSize',fontSize)
xlabel('$q$','Interpreter','latex','FontSize',fontSize);
ylabel('$Counts$','Interpreter','latex','FontSize',fontSize);
xlim(axesHistogram1,[0 histHeader1.BinLimits(2)]);
set(axesHistogram1,'FontName','arial','FontSize',24,'LineWidth',1);

% %% PLOT THE UNSCALED HISTOGRAM

% FIT THE LASING UNSCALED DATA

lasngPd = fitdist(DATA,'Kernel','Kernel','normal');
%lasngPd = fitdist(DATA,'Stable');
% 
histFigure2 = figure('InvertHardcopy','off','Color',[1 1 1]);

% AXES
axesHistogram2 = axes('Parent',histFigure2); hold(axesHistogram2,'on');

histHeader2 = histogram(DATA,'Parent',axesHistogram2,'LineWidth',1,...
            'BinMethod','auto','Normalization','PDF');

xFit2 = 0:histHeader2.BinWidth:histHeader2.BinLimits(2);
yFit2 = pdf(lasngPd,xFit2);
%mixxFit2 = xFit2;
%mixyfit2 = yFit2;
realxFit2 = xFit2;
realyfit2 = yFit2;
line(xFit2,yFit2,'LineWidth',2,'Color',[1.00,0.41,0.16]);
title(['FLUENCE ' strPOWER '\muJ/cm^2']);
ylabel({'PDF'});
xlabel({'unscaled q (a.u.)'});
xlim(axesHistogram2,[0 histHeader2.BinLimits(2)]);
%ylim(axes1,[0 3]);
box(axesHistogram2,'on');
set(axesHistogram2,'FontName','arial','FontSize',24,'LineWidth',1);

FIG = [dataPath dataFolder '/rawLasngHist_pow' strPOWER]; 

savefig(gcf,[FIG '.fig']);
print(gcf,[FIG '.eps'],'-deps','-painters')
print(gcf,[FIG '.png'],'-dpng','-painters')

%% CREATE A SYNTHETIC DATASET
% 
% fakeLasngQ = random(lasngPd,length(DATA));
% 
% FILENAME = [dataPath dataFolder '/fakeLasngQ_' num2str(Data.power(p)) '[uJcm2].csv'];
% writematrix(fakeLasngQ,FILENAME,'Delimiter',',');
% 
% PLOT THE UNSCALED FAKE DATA
% 
% 
% histFigure3 = figure('InvertHardcopy','off','Color',[1 1 1]);
% 
% AXES
% axesHistogram3 = axes('Parent',histFigure3); hold(axesHistogram3,'on');
% title(['fake q data at ' strPOWER '\muJ/cm^2']);
% 
% histHeader3 = histogram(fakeLasngQ,'Parent',axesHistogram3,'LineWidth',1,...
%             'BinMethod','auto','Normalization','PDF');
% 
% xFit3 = 0:histHeader3.BinWidth:histHeader3.BinLimits(2);
% yFit3 = pdf(lasngPd,xFit3);
% line(xFit3,yFit3,'LineWidth',2,'Color',[1.00,0.41,0.16]);
% 
% ylabel({'PDF'});
% xlabel({'unscaled q (a.u.)'});
% xlim(axesHistogram3,[0 histHeader3.BinLimits(2)]);
% ylim(axes1,[0 3]);
% box(axesHistogram3,'on');
% set(axesHistogram3,'FontName','arial','FontSize',24,'LineWidth',1);
% 
% FIG = [dataPath dataFolder '/fakeLasngHist_pow' strPOWER]; 
% 
% savefig(gcf,[FIG '.fig']);
% print(gcf,[FIG '.eps'],'-deps','-painters')
% print(gcf,[FIG '.png'],'-dpng','-painters')
% 
% 
% hold(axesHistogram3,'off');


end

