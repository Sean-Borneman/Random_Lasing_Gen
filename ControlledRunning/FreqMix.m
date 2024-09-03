%% Clean UP
close all;
clear all;
clc;
rng("default")
%% Load Data
load('SpectralStatData.mat');
POWER = 8;

numberOfReplica = size(Data.singleSpec(1).data,2);
numberOfWavelen = size(Data.singleSpec(1).data,1);
numberOfDataPow = length(Data.power); 
newPulseCollection = zeros(1340, 100);

dataPath = 'C:/Users/Sean/Documents/Programing/Python/PURDUE/ControlledRunning/GeneratedDataGradual';
%% Iterate over each power

for p =7:7% 1:numberOfDataPow %1:1 7:7%
    pulseCollection = Data.singleSpec(p).data;
    numPOWER = round(Data.power(p)*1000);
    %% Reshape collection
    rows = size(pulseCollection);
    rows = rows(1);
    
    % Initialize a cell array to hold the 1x3 lists
    array_of_lists = cell(1, rows);

    % Extract each row and store it in the cell array
    for i = 1:rows
        array_of_lists{i} = pulseCollection(i, :);
    end

    new_pulse =  zeros(1340, 100);%[];
    last_Selection = -1;
    for pulses = 1:100
        new_shot = zeros(1340, 1);
        for wavelength = 1:numberOfWavelen
            if last_Selection ~= -1
                selectionList = filtList(array_of_lists{wavelength}, last_Selection-10, last_Selection+10);
                

                if length(selectionList) > 0
                    selection = randsample(selectionList,1);
                else
                    selection = randsample(pulseCollection(wavelength, :),1);
                end
            else
                selection = randsample(pulseCollection(wavelength, :),1);
            end
            index = find(array_of_lists{wavelength} == selection, 1);
            array_of_lists{wavelength}(index) = [];
            
            new_shot(wavelength) = selection;
            
            last_Selection = selection;
%             indx = find(pulseCollection == selection, 1); 
%             pulseCollection(indx) = [];
        end
        new_pulse(:, pulses) = new_shot;
        
    end
    
    figure;
    imagesc(Data.singleSpec(p).data);
    
    figure 
    imagesc(new_pulse);
    
%     [genFitX,genFitY] = getFits(new_pulse);
%     [realFitX,realFitY] = getFits(Data.singleSpec(p).data);
%     disp("Loss:");
%     disp(getLoss(genFitX, genFitY, realFitX, realFitY));
%%Lasing Loss
%     [genFitX,genFitY] = getLasingFits(new_pulse);
%     [realFitX,realFitY] = getLasingFits(Data.singleSpec(p).data);
%     lasing_loss = getLoss(genFitX, genFitY, realFitX, realFitY, "Lasing");
%     errors = [errors, lasing_loss];
    
    %Spont Loss
    [genFitX,genFitY] = getSpontFits(new_pulse);
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

function filteredList = filtList(IList, lower_bound, upper_bound)
    filteredList = [];
    for item = 1:length(IList)
        if IList(item) < upper_bound && IList(item) > lower_bound
            filteredList = [filteredList, IList(item)];
        end
    end
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

function loss = getLoss(genFitX, genFitY, realFitX, realFitY, titleN)
    figure
    hold on
    line(genFitX,genFitY,'LineWidth',2,'Color',[1.00,0.41,0.16]);
    line(realFitX,realFitY,'LineWidth',2,'Color','blue', 'LineStyle','--');
    title("Compared Orderr Parameters for" + titleN + "reigem");
    legend('Generated Data', 'Real Data')
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
% 
% %% Clean UP
% close all;
% clear all;
% clc;
% rng("default")
% %% Load Data
% load('SpectralStatData.mat');
% POWER = 8;
% 
% numberOfReplica = size(Data.singleSpec(1).data,2);
% numberOfWavelen = size(Data.singleSpec(1).data,1);
% numberOfDataPow = length(Data.power); 
% newPulseCollection = zeros(1340, 100);
% 
% dataPath = 'C:/Users/Sean/Documents/Programing/Python/PURDUE/ControlledRunning/GeneratedDataGradual';
% %% Iterate over each power
% 
% for p = 7:7%1:numberOfDataPow
%     pulseCollection = Data.singleSpec(p).data;
%     %% Reshape collection
%     rows = size(pulseCollection);
%     rows = rows(1);
%     
%     % Initialize a cell array to hold the 1x3 lists
%     array_of_lists = cell(1, rows);
% 
%     % Extract each row and store it in the cell array
%     for i = 1:rows
%         array_of_lists{i} = pulseCollection(i, :);
%     end
% 
%     new_pulse =  zeros(1340, 100);%[];
%     last_Selection = -1;
%     for pulses = 1:100
%         new_shot = zeros(1340, 1);
%         for wavelength = 1:numberOfWavelen
%             if last_Selection ~= -1
%                 selectionList = filtList(pulseCollection(wavelength, :), last_Selection-30, last_Selection+30);
%                 if length(selectionList) > 0
%                     selection = randsample(selectionList,1);
%                 else
%                     selection = randsample(pulseCollection(wavelength, :),1);
%                 end
%             else
%                 selection = randsample(pulseCollection(wavelength, :),1);
%             end
%             new_shot(wavelength) = selection;
%             
%             last_Selection = selection;
% %             indx = find(pulseCollection == selection, 1); 
% %             pulseCollection(indx) = [];
%         end
%         new_pulse(:, pulses) = new_shot;
%         
%     end
%     
%     figure;
%     imagesc(Data.singleSpec(p).data);
%     
%     figure 
%     imagesc(new_pulse);
%     
% %     [genFitX,genFitY] = getFits(new_pulse);
% %     [realFitX,realFitY] = getFits(Data.singleSpec(p).data);
% %     disp("Loss:");
% %     disp(getLoss(genFitX, genFitY, realFitX, realFitY));
% end
% 
% function filteredList = filtList(IList, lower_bound, upper_bound)
%     filteredList = [];
%     for item = 1:length(IList)
%         if IList(item) < upper_bound && IList(item) > lower_bound
%             filteredList = [filteredList, IList(item)];
%         end
%     end
% end
% 
% function [xFited, yFited] = getFits(pulses)
%     load('SpectralStatData.mat');
%     % 
%     %freqMix =table2array(readtable('freqMIX238.csv'));% table2array(readtable('avrgwDist.csv'));
%     Data.singleSpec(1).data = pulses
%     %% CONSTANTS
%     close all  
%     binsNumber = 30;
%     fontSize   = 20;
%     wavelength = Data.Wavel;
%     lasngIndex = (find(wavelength<=550,1,'last'):find(wavelength>=557,1)); % 7 nm centered around lasing wavel
%     spontIndex = (find(wavelength<=536,1,'last'):find(wavelength>=543,1)); % 7 nm centered around 
%     powerIndex = (1:length(Data.power));
%     lasngCentr = [556,555,554.5,554,553.5,553,552.5,552,551.5,551,551,551,551,551];
%     %% MAIN LOOP - SCAN ALL THE 
%     numberOfReplica = size(Data.singleSpec(1).data,2);
%     numberOfWavelen = size(Data.singleSpec(1).data,1);
%     numberOfDataPow = length(Data.power); 
%     sizeQ = numberOfDataPow*numberOfDataPow-numberOfDataPow;
%     flatSpontQ     = zeros(sizeQ,1);
%     flatLasngQ     = zeros(sizeQ,1);
%     flatNormSpontQ = zeros(sizeQ,1);
%     flatNormLasngQ = zeros(sizeQ,1);
%     for p = 1:1%numberOfDataPow
% 
%     close all;
% 
%     numPOWER = round(Data.power(p)*1000);
%     strPOWER = num2str(numPOWER);
% 
%     dataFolder = ['Power_' strPOWER '[uJcm2]']; 
% 
%     mkdir(dataFolder);
% 
%     %% STEP 1: Take the 7-nm band centered around the lasing wavelength
%     lasngIndex = (find(wavelength<=(lasngCentr(p)-3.5),1,'last'):find(wavelength>=(lasngCentr(p)+3.5),1));
% 
%     %% STEP 2: GET AND STORE SPECTROGAMM MATRIX Ip FOR A GIVEN FLUENCE P  
%     Ip      = Data.singleSpec(p).data;
%     spontIp = Ip(spontIndex,:);
%     lasngIp = Ip(lasngIndex,:);
% 
%     %% STEP 3: AVERAGE ALL REPLICAS IN Ip TO GET A MEAN VECTOR, 
%     %% REPLICATED INTO A MATRIX Ibar AND STORE IT IN A FILE
%     spontIbar = repmat(mean(spontIp,2),[1,size(spontIp,2)]);
%     lasngIbar = repmat(mean(lasngIp,2),[1,size(lasngIp,2)]);
% 
%     %% STEP 4: GET THE DIFFERENCE MATRICES \Delta FOR A GIVEN FLUENCE P 
%     spontDelta     = spontIp - spontIbar;
%     lasngDelta     = lasngIp - lasngIbar;
% 
% 
%     %% STEP 5: NORMALIZE THE DIFFERENCE MATRICES \Delta 
%     wavelengthNumber     = size(spontDelta,1);
%     absSpontDelta        = repmat(sqrt(sum(spontDelta.^2,1)),[wavelengthNumber,1]);
%     normalizedSpontDelta = spontDelta./absSpontDelta; 
% 
%     wavelengthNumber     = size(lasngDelta,1);
%     absLasngDelta        = repmat(sqrt(sum(lasngDelta.^2,1)),[wavelengthNumber,1]);
%     normalizedLasngDelta = lasngDelta./absLasngDelta;
% 
%     spontQ = abs(spontDelta'*spontDelta);
%     lasngQ = abs(lasngDelta'*lasngDelta);
% 
%     normalizedSpontQ = abs(normalizedSpontDelta'*normalizedSpontDelta);
%     normalizedLasngQ = abs(normalizedLasngDelta'*normalizedLasngDelta);
% 
% 
%     %% IGNORE THE DATA AT THE DIAGONAL
%     %  AND FLATTEN THE ARRAY
%     a  = numberOfReplica;
%     a1 = a+1;
%     a2 = a*a;
% 
%     flatSpontQ = reshape(spontQ, a2,1);
%     flatLasngQ = reshape(lasngQ, a2,1);
% 
%     flatNormSpontQ = reshape(normalizedSpontQ,a2,1);
%     flatNormLasngQ = reshape(normalizedLasngQ,a2,1);
% 
%     % flatSpontQ(1:a1:a2)=[];
%     % flatLasngQ(1:a1:a2)=[];
%     % 
%     % flatNormSpontQ(1:a1:a2)=[];
%     % flatNormLasngQ(1:a1:a2)=[];
% 
%     Ibar = repmat(mean(Ip,2),[1,2]);
% 
%     x = [0, 101, 101, 0, 0];
%     y = [wavelength(lasngIndex(1)),wavelength(lasngIndex(1)),wavelength(lasngIndex(end)),wavelength(lasngIndex(end)),wavelength(lasngIndex(1))];
%     %% RAW DELTAS (w/o normalization)
% 
%     matrixA = abs(spontDelta);
%     matrixB = abs(lasngDelta);
% 
%     %% NORMALIZED DELTAS
% 
%     matrixA = abs(normalizedSpontDelta);
%     matrixB = abs(normalizedLasngDelta);
%     %% RAW Q (w/o normalization)
% 
%     matrixA = abs(spontQ);
%     matrixB = abs(lasngQ);
% 
%     
%     %% NORMALIZED Q
% 
%     matrixA = abs(normalizedSpontQ);
%     matrixB = abs(normalizedLasngQ);
% 
%     %% Raw Histograms
%     close all
% 
%     matrixA = abs(flatSpontQ);
%     matrixB = abs(flatLasngQ);
% 
% 
%     
%     %% NORMALIZED HISTOGRAM
% 
%     matrixA = abs(flatNormSpontQ);
%     matrixB = abs(flatNormLasngQ);
% 
%     
% 
%     %% other Histogram
%     DATA = [flatLasngQ; -flatLasngQ];
% 
% 
%     histFigure1 = figure('InvertHardcopy','off','Color',[1 1 1]);
%     % AXES
%     axesHistogram1 = axes('Parent',histFigure1); hold(axesHistogram1,'on');
%     histHeader1 = histogram(DATA,'EdgeColor','k','EdgeAlpha',0.3,'FaceColor','m');
% 
%     % FIT THE LASING UNSCALED DATA
% 
%     lasngPd = fitdist(DATA,'Kernel','Kernel','normal');
%     %lasngPd = fitdist(DATA,'Stable');
%     % 
%     histFigure2 = figure('InvertHardcopy','off','Color',[1 1 1]);
% 
%     % AXES
%     axesHistogram2 = axes('Parent',histFigure2); hold(axesHistogram2,'on');
% 
%     histHeader2 = histogram(DATA,'Parent',axesHistogram2,'LineWidth',1,...
%                 'BinMethod','auto','Normalization','PDF');
% 
%     xFited = 0:histHeader2.BinWidth:histHeader2.BinLimits(2);
%     yFited = pdf(lasngPd,xFited);
%     %mixxFit2 = xFit2;
%     %mixyfit2 = yFit2;
% %     realxFit2 = xFit2;
% %     realyfit2 = yFit2;
%     end
% end
% 
% function error = getAllowedVarience(newPulseCollection, f, i)
%     x = 1:30;
%     y = newPulseCollection(f-30:f-1, i);
%     y = reshape(y,1,30);
%     fit = polyfit(x,y , 1);
%     y1 = polyval(fit, x);
% %     figure
% %     plot(x,y,'o')
% %     title("Best Fit")
% %     xlabel("Wavelength");
% %     ylabel("Amplitude");
% %     hold on
% %     plot(x,y1)
% %     hold off
%     error = abs(y1(1)-y1(length(y1)))/length(y1);%+0.028*y(1);  
% end
% 
% function loss = getLoss(genFitX, genFitY, realFitX, realFitY)
%     line(genFitX,genFitY,'LineWidth',2,'Color',[1.00,0.41,0.16]);
%     line(realFitX,realFitY,'LineWidth',2,'Color',[1.00,0.41,0.16]);
% 
%     d1 = realFitY;
%     d1_len = numel(realFitY);
%     d2 = genFitY;
%     d2_len = numel(genFitY);
%     try
%         d1 = interp1(1:(d2_len / d1_len):d2_len,d1,1:d2_len,'linear','extrap');
%     catch
%         try
%         d1 = interp1(1:(d2_len / d1_len):d2_len+(d2_len / d1_len),d1,1:d2_len,'linear','extrap');
%         catch
%             d1 = interp1(1:(d2_len / d1_len):d2_len+2*(d2_len / d1_len),d1,1:d2_len,'linear','extrap');
%         end
%     end
%     %d1 = interp1(1:(d2_len / d1_len):d2_len,d1,1:d2_len,'linear','extrap');
%     %line(mixxFit2,d1)
%     RMSE = sqrt(mean(((genFitY - d1) .^ 2)));
%     loss = RMSE;
% end
% 
% 
% 
% 
% 
