%% Clean UP
close all;
clear all;
clc;
rng("default")
%% Load Data
load('SpectralStatData.mat');
POWER = 8;
errors = [];
numberOfReplica = size(Data.singleSpec(1).data,2);
numberOfWavelen = size(Data.singleSpec(1).data,1);
numberOfDataPow = length(Data.power); 
newPulseCollection = zeros(1340, 100);
dataPath = 'C:/Users/Sean/Documents/Programing/Python/PURDUE/ControlledRunning/GeneratedDataGradual';
%% Iterate over each power

for p = 1:numberOfDataPow %7:7
    %% Graph each shot Collection
    numPOWER = round(Data.power(p)*1000);
    disp("POWER______:"+numPOWER);
    strPOWER = num2str(numPOWER);
    dataFolder = [dataPath '/Power_' strPOWER '[uJcm2]']; 
    mkdir(dataFolder);
    
%     figure;
%     imagesc(Data.singleSpec(p).data);
%     title("Data from pulse"+ strPOWER+"[ujcm2]");
%     saveas(gcf,strcat(num2str(dataFolder), '/InitialData.png'))
%     
    pulseCollection = Data.singleSpec(p).data;
    frequencyHist = pulseCollection(925, :) ;
    
%     figure
%     histogram(frequencyHist(:))
%     title("Amplitude Histogram at " +Data.Wavel(925)+"[nm]");
    for f = 1:length(Data.Wavel)
        pulseCollection = Data.singleSpec(p).data;
%         figure
%         plot(pulseCollection(:, f));
%         title("Wavelength Amplitude");
    

        %fit theh first frequency to a normal distribution, must flip bc need a collumn not row vector
        frequencyHist = pulseCollection(f, :) ;
        pd = fitdist(frequencyHist(:),'Kernel','Kernel','epanechnikov');%fitdist(frequencyHist(:),'Gamma');
    
%         figure
%         histfit(frequencyHist(:),10, 'Kernel');
%         %show the amplitude distribution
% %         figure
% %         histogram(frequencyHist(:))
%          title("Amplitude Histogram at " +Data.Wavel(925)+"[nm]");
%          xlabel("Amplitude");
%          ylabel("counts");
%         hold on
%         %fit the distribution to a normal PDF
%         x_values = 122:.1:140;
%         y = pdf(pd,x_values);
%         
% %         figure;
%         plot(x_values,y);
%         hold off
%         title("PDF for frequency");
        %% Generate 100 amplitudes for this wavelength
        %samples a random number from the PDF
        for i = 1:100
            newPoint = random(pd);
            if f > 31
                ticks = 0;
                
                AllowedVarience = getAllowedVarience(newPulseCollection, f, i);
                flatError = 50;
                while abs(newPulseCollection(f-1, i) - newPoint) > AllowedVarience+flatError%10%+10*newPoint^0.35
                    newPoint = random(pd);
                    ticks=ticks+1;
                    if ticks > 10
                        
                        flatError = flatError + 5;
                    end
                end
                newPulseCollection(f, i) = newPoint;
            else
                newPulseCollection(f, i) = random(pd);
                
            end
        
            
        end
    %disp(f);
    end
%     figure;
%     imagesc(newPulseCollection);
%     title("Generated Data from pulse"+ strPOWER+"[ujcm2]");
%     saveas(gcf,strcat(num2str(dataFolder), '/GeneratedData.png'))
%     writematrix(newPulseCollection,strcat(dataFolder, "/GeneratedShots.csv"), 'Delimiter',',')
    %FILENAME = "RawData"+ round(Data.power(POWER)*1000)+ "[uJcm2].csv";%
    %writematrix(Data.singleSpec(POWER).data,FILENAME,'Delimiter',',');
%     figure
%     plot(newPulseCollection(:, 1));
%     title("Wavelength Amplitude atGen");
%     saveas(gcf,strcat(num2str(dataFolder), '/SamplePulse.png'))
%     figure
%     plot(Data.singleSpec(p).data(:, 1));
%     title("Wavelength Amplitude at Pulse");
%     saveas(gcf,strcat(num2str(dataFolder), '/SamplePulse.png'))
%     
%     figure;
%     imagesc(Data.singleSpec(p).data);
%     title("Data from pulse"+ strPOWER+"[ujcm2]");
%     
    
    [genFitX,genFitY] = getFits(newPulseCollection);
    [realFitX,realFitY] = getFits(Data.singleSpec(p).data);
    disp("Loss:" + flatError);
    final_loss = getLoss(genFitX, genFitY, realFitX, realFitY);
    disp(final_loss);
    errors = [errors, final_loss];
end
disp(errors);

function [xFited, yFited] = getFits(pulses)
    load('SpectralStatData.mat');
    % 
    %freqMix =table2array(readtable('freqMIX238.csv'));% table2array(readtable('avrgwDist.csv'));
    Data.singleSpec(1).data = pulses
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

    %% STEP 3: AVERAGE ALL REPLICAS IN Ip TO GET A MEAN VECTOR, 
    %% REPLICATED INTO A MATRIX Ibar AND STORE IT IN A FILE
    spontIbar = repmat(mean(spontIp,2),[1,size(spontIp,2)]);
    lasngIbar = repmat(mean(lasngIp,2),[1,size(lasngIp,2)]);

    %% STEP 4: GET THE DIFFERENCE MATRICES \Delta FOR A GIVEN FLUENCE P 
    spontDelta     = spontIp - spontIbar;
    lasngDelta     = lasngIp - lasngIbar;


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

    Ibar = repmat(mean(Ip,2),[1,2]);

    x = [0, 101, 101, 0, 0];
    y = [wavelength(lasngIndex(1)),wavelength(lasngIndex(1)),wavelength(lasngIndex(end)),wavelength(lasngIndex(end)),wavelength(lasngIndex(1))];
    %% RAW DELTAS (w/o normalization)

    matrixA = abs(spontDelta);
    matrixB = abs(lasngDelta);

    %% NORMALIZED DELTAS

    matrixA = abs(normalizedSpontDelta);
    matrixB = abs(normalizedLasngDelta);
    %% RAW Q (w/o normalization)

    matrixA = abs(spontQ);
    matrixB = abs(lasngQ);

    
    %% NORMALIZED Q

    matrixA = abs(normalizedSpontQ);
    matrixB = abs(normalizedLasngQ);

    %% Raw Histograms
    close all

    matrixA = abs(flatSpontQ);
    matrixB = abs(flatLasngQ);


    
    %% NORMALIZED HISTOGRAM

    matrixA = abs(flatNormSpontQ);
    matrixB = abs(flatNormLasngQ);

    

    %% other Histogram
    DATA = [flatLasngQ; -flatLasngQ];


    histFigure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
    % AXES
    axesHistogram1 = axes('Parent',histFigure1); hold(axesHistogram1,'on');
    histHeader1 = histogram(DATA,'EdgeColor','k','EdgeAlpha',0.3,'FaceColor','m');

    % FIT THE LASING UNSCALED DATA

    lasngPd = fitdist(DATA,'Kernel','Kernel','normal');
    %lasngPd = fitdist(DATA,'Stable');
    % 
    histFigure2 = figure('InvertHardcopy','off','Color',[1 1 1]);

    % AXES
    axesHistogram2 = axes('Parent',histFigure2); hold(axesHistogram2,'on');

    histHeader2 = histogram(DATA,'Parent',axesHistogram2,'LineWidth',1,...
                'BinMethod','auto','Normalization','PDF');

    xFited = 0:histHeader2.BinWidth:histHeader2.BinLimits(2);
    yFited = pdf(lasngPd,xFited);
    %mixxFit2 = xFit2;
    %mixyfit2 = yFit2;
%     realxFit2 = xFit2;
%     realyfit2 = yFit2;
    end
end

function error = getAllowedVarience(newPulseCollection, f, i)
    x = 1:30;
    y = newPulseCollection(f-30:f-1, i);
    y = reshape(y,1,30);
    fit = polyfit(x,y , 1);
    y1 = polyval(fit, x);
%     figure
%     plot(x,y,'o')
%     title("Best Fit")
%     xlabel("Wavelength");
%     ylabel("Amplitude");
%     hold on
%     plot(x,y1)
%     hold off
    error = abs(y1(1)-y1(length(y1)))/length(y1);%+0.028*y(1);  
end

function loss = getLoss(genFitX, genFitY, realFitX, realFitY)
    line(genFitX,genFitY,'LineWidth',2,'Color',[1.00,0.41,0.16]);
    line(realFitX,realFitY,'LineWidth',2,'Color',[1.00,0.41,0.16]);

    d1 = realFitY;
    d1_len = numel(realFitY);
    d2 = genFitY;
    d2_len = numel(genFitY);
    try
        d1 = interp1(1:(d2_len / d1_len):d2_len,d1,1:d2_len,'linear','extrap');
    catch
        try
        d1 = interp1(1:(d2_len / d1_len):d2_len+(d2_len / d1_len),d1,1:d2_len,'linear','extrap');
        catch
            d1 = interp1(1:(d2_len / d1_len):d2_len+2*(d2_len / d1_len),d1,1:d2_len,'linear','extrap');
        end
    end
    %d1 = interp1(1:(d2_len / d1_len):d2_len,d1,1:d2_len,'linear','extrap');
    %line(mixxFit2,d1)
    RMSE = sqrt(mean(((genFitY - d1) .^ 2)));
    loss = RMSE;
end



