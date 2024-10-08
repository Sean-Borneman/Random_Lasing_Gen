%% Clean UP
close all;
clear all;
clc;
rng("default")
warning('off','all')
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

for p = 7:7%1:numberOfDataPow %7:7 
    %% Graph each shot Collection
    numPOWER = round(Data.power(p)*1000);
    disp("POWER______:"+numPOWER);
    strPOWER = num2str(numPOWER);
%     dataFolder = [dataPath '/Power_' strPOWER '[uJcm2]']; 
%     mkdir(dataFolder);
    pulseCollection = Data.singleSpec(p).data;
    averagePulse = mean(pulseCollection, 2);
    diffCollection = pulseCollection - averagePulse;
    
     figure;
     imagesc(diffCollection);
     title("Diff Data from pulse"+ strPOWER+"[ujcm2]");
    
%     known bug where when counting backwards the slope based error is 0
    for f =930:-1:1% 1:length(Data.Wavel) 
%         fit theh first frequency to a normal distribution, must flip bc need a collumn not row vector
        frequencyHist = diffCollection(f, :) ;
%         Fit to stable distribution, Levi distribution is a special case
%         where alpha = 0.5, beta = 1
        pd = fitdist(frequencyHist(:),'Stable');%fitdist(frequencyHist(:),'Gamma');
        if (f == 930)
            %pd.gam = pd.gam +10;
            data = diffCollection(f, :);
%             Generate x values for plotting the PDF
            x_values = linspace(min(data), max(data), 1000);

%             Calculate PDF values of the fitted normal distribution
            pdf_values = pdf(pd, x_values);

%             Plot the PDF
            figure;
            hold on
            histogram(data, 25,'Normalization', 'pdf');
            plot(x_values, pdf_values, 'r-', 'LineWidth', 2);
            xlabel('Data');
            ylabel('Probability Density');
            title('Normal Distribution PDF');
            legend('Fitted Normal PDF');
            hold off
        end
        % Generate 100 amplitudes for this wavelength
%         samples a random number from the PDF
        for i = 1:100
            newPoint = random(pd);
            if f > 31
                ticks = 0;
                
                AllowedVarience = getAllowedVarienceBackward(newPulseCollection, f, i);
                flatError = 0;
                if f == 930
                    flatError = 5000;
                end
                if f <= 825
                    flatError = 500;
                end
                
                while abs(newPulseCollection(f+1, i) - newPoint) > AllowedVarience+flatError%10%+10*newPoint^0.35
                    newPoint = random(pd);
                    ticks=ticks+1;
                    if ticks > 10
                        ticks = 0;
                        flatError = flatError + 0.1; % was 0.1 %0.3 gest me close
                    end
                end
                newPulseCollection(f, i) = newPoint;
            else
                newPulseCollection(f, i) = random(pd);
                
            end 
        end
        disp(f);
    end
    for f =931:length(Data.Wavel) 
%         fit theh first frequency to a normal distribution, must flip bc need a collumn not row vector
        frequencyHist = diffCollection(f, :) ;
%         Fit tro stable distribution, Levi distribution is a special case
%         where alpha = 0.5, beta = 1
        pd = fitdist(frequencyHist(:),'Stable');%fitdist(frequencyHist(:),'Gamma');
        if (f == 930)
            data = diffCollection(f, :);
%             Generate x values for plotting the PDF
            x_values = linspace(min(data), max(data), 1000);

%             Calculate PDF values of the fitted normal distribution
            pdf_values = pdf(pd, x_values);

%             Plot the PDF
            figure;
            hold on
            histogram(data, 25,'Normalization', 'pdf');
            plot(x_values, pdf_values, 'r-', 'LineWidth', 2);
            xlabel('Data');
            ylabel('Probability Density');
            title('Normal Distribution PDF');
            legend('Fitted Normal PDF');
        end
        % Generate 100 amplitudes for this wavelength
%         samples a rando m number from the PDF
        for i = 1:100
            newPoint = random(pd);
            if f > 31
                ticks = 0;
                
                AllowedVarience = getAllowedVarience(newPulseCollection, f, i);
                flatError = 0;
                if f == 930
                    flatError = 5000;
                end
                
                while abs(newPulseCollection(f-1, i) - newPoint) > AllowedVarience+flatError%10%+10*newPoint^0.35
                    newPoint = random(pd);
                    ticks=ticks+1;
                    if ticks > 10
                        ticks = 0;
                        flatError = flatError + 0.3;%was 0.1
                    end
                end
                newPulseCollection(f, i) = newPoint;
            else
                newPulseCollection(f, i) = random(pd);
                
            end 
        end
        disp(f);
    end
    newPulseCollection = newPulseCollection + averagePulse + wgn(1340, 100, 23);%was 30
%     newPulseCollection = averagePulse + wgn(1340, 100, 23);
    clims = [min(pulseCollection, [], 'all') max(pulseCollection, [], 'all')];
    %set Color bounds
    figure;
    hold on
    imagesc(newPulseCollection, clims);
    title("Gen Data from pulse"+ strPOWER+"[ujcm2]");
    hold off
    %%Lasing Loss
%     [genFitX,genFitY] = getLasingFits(newPulseCollection);
%     [realFitX,realFitY] = getLasingFits(Data.singleSpec(p).data);
%     lasing_loss = getLoss(genFitX, genFitY, realFitX, realFitY, "Lasing");
%     errors = [errors, lasing_loss];
    
    %%Spont Loss
    [genFitX,genFitY] = getSpontFits(newPulseCollection);
    [realFitX,realFitY] = getSpontFits(Data.singleSpec(p).data);
    spont_loss = getLoss(genFitX, genFitY, realFitX, realFitY, "Spont");
    disp(spont_loss);
    disp("POWER______:"+numPOWER);
%     disp("Lasing Loss:");
%     disp(lasing_loss);
    disp("Spont Loss:");
    disp(spont_loss);
%     disp("Total Error:");
%     disp(lasing_loss + spont_loss);
end

function [xFited, yFited] = getLasingFits(pulses)
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
function error = getAllowedVarienceBackward(newPulseCollection, f, i)
    x = 1:30;
    y = newPulseCollection(f+1:f+30, i);
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
function loss = getLoss(genFitX, genFitY, realFitX, realFitY, titleN)
    figure
    hold on 
    line(genFitX,genFitY,'LineWidth',2,'Color',[1.00,0.41,0.16]);
    line(realFitX,realFitY,'LineWidth',2,'Color','blue', 'LineStyle','--');
    title("Compared Orderr Parameters for" + titleN + "reigem");
    legend('on');
    hold off
    
    d1 = realFitY;
    d1_len = numel(realFitY);
    d2 = genFitY;
    d2_len = numel(genFitY);
    
    d1 = interp1(1:d1_len, d1, d2_len,'linear', 'extrap'); %doesnt work becuase d2 can query points outside of d1's function boundry
%     try
%         d1 = interp1(1:(d2_len / d1_len):d2_len,d1,1:d2_len,'linear','extrap');
%     catch
%         try
%         d1 = interp1(1:(d2_len / d1_len):d2_len+(d2_len / d1_len),d1,1:d2_len,'linear','extrap');
%         catch
%             d1 = interp1(1:(d2_len / d1_len):d2_len+2*(d2_len / d1_len),d1,1:d2_len,'linear','extrap');
%         end
%     end
    %d1 = interp1(1:(d2_len / d1_len):d2_len,d1,1:d2_len,'linear','extrap');
    %line(mixxFit2,d1)
    RMSE = sqrt(mean(((genFitY - d1) .^ 2)));
    loss = RMSE;
end

function [xFited, yFited] = getSpontFits(pulses)
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
    DATA = [flatSpontQ; -flatSpontQ];


    histFigure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
    % AXES
    axesHistogram1 = axes('Parent',histFigure1); hold(axesHistogram1,'on');
    histHeader1 = histogram(DATA,'EdgeColor','k','EdgeAlpha',0.3,'FaceColor','m');

    % FIT THE LASING UNSCALED DATA

    spontPd = fitdist(DATA,'Kernel','Kernel','normal');
    %lasngPd = fitdist(DATA,'Stable');
    % 
    histFigure2 = figure('InvertHardcopy','off','Color',[1 1 1]);

    % AXES
    axesHistogram2 = axes('Parent',histFigure2); hold(axesHistogram2,'on');

    histHeader2 = histogram(DATA,'Parent',axesHistogram2,'LineWidth',1,...
                'BinMethod','auto','Normalization','PDF');

    xFited = 0:histHeader2.BinWidth:histHeader2.BinLimits(2);
    yFited = pdf(spontPd,xFited);
    %mixxFit2 = xFit2;
    %mixyfit2 = yFit2;
%     realxFit2 = xFit2;
%     realyfit2 = yFit2;
    end
end

