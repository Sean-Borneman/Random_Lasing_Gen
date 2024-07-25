close all 
clear all
clc
%Fit probability Distribution to data

load examgrades
x = grades(:, 1);
x_pdf = [1:0.1:100];


%---NORMAL------
disp("Normal")
a = 40;
b = 10;
dist = sqrt(a^2 + b^2);
pd_true = makedist('Normal','mu',a,'sigma',b);

%rng('default') % For reproducibility
tr = 0;
error_list = [20:1:1000];
rmse_list = [20:1:1000];
for sn = 20:1000
    y_true = reshape([1:1:sn], [], 1);
    for i = 1:sn
        y_true(i) = random(pd_true);
    end

    pd = fitdist(y_true, 'Normal');

    error = sqrt((pd_true.mu-pd.mu)^2+(pd_true.sigma-pd.sigma)^2);
    tr = tr +error;
    error_list(i-19) = error/dist;
    
    y = pdf(pd, x_pdf);
    y2 = pdf(pd_true, x_pdf);
    ss = 0;
    for w=1:length(y)
        ss = ss + (y2(w)-y(w))^2;
    end
    rmse = sqrt(ss/length(y));
    rmse_list(i-19) = rmse;
end
figure
plot([20:1000],error_list)
title("Normal Distribution: Loss across Samples")
xlabel('Number of Samples') 
ylabel('Percent Loss') 
disp("average error:" + tr/981)
disp("percent Error:" + (tr/981)/dist)

figure
plot([20:1000],rmse_list)
title("Normal Distribution: Loss across Samples")
xlabel('Number of Samples') 
ylabel('RMSE Loss') 


%%--- Gamma------%%

disp("Gamma")
a = 40;
b = 10;
dist = sqrt(a^2 + b^2);
pd_true = makedist('Gamma','a',a,'b',b);

rng('default') % For reproducibility
tr = 0;
error_list = [20:1:1000];
rmse_list = [20:1:1000];


for sn = 20:1000
    y_true = reshape([1:1:sn], [], 1);
    for i = 1:sn
        y_true(i) = random(pd_true);
    end

    pd = fitdist(y_true+awgn(y_true,1)*0.01, 'Gamma');

    error = sqrt((pd_true.a-pd.a)^2+(pd_true.b-pd.b)^2);
    tr = tr +error;
    error_list(i-19) = error/dist;
    %y = pdf(pd, x_pdf);
    y = pdf(pd, x_pdf);
    y2 = pdf(pd_true, x_pdf);
    ss = 0;
    for w=1:length(y)
        ss = ss + (y2(w)-y(w))^2;
    end
    rmse = sqrt(ss/length(y));
    rmse_list(i-19) = rmse;
end
figure
plot([20:1000],error_list)
title("Gamma Distribution: Loss across Samples")
xlabel('Number of Samples') 
ylabel('Percent Loss') 
disp("average error:" + tr/981)
disp("percent Error:" + (tr/981)/dist)

figure
plot([20:1000],rmse_list)
title("Gamma Distribution: Loss across Samples")
xlabel('Number of Samples') 
ylabel('RMSE Loss') 
%%--- EXP------%%

disp("Exponential")
a = 40;
dist = sqrt(a^2);
pd_true = makedist('Exponential','mu',a);

rng('default') % For reproducibility
tr = 0;
error_list = [20:1:1000];
rmse_list = [20:1:1000];


for sn = 20:1000
    y_true = reshape([1:1:sn], [], 1);
    for i = 1:sn
        y_true(i) = random(pd_true);
    end

    pd = fitdist(abs(y_true+awgn(y_true,1)*0.01), 'Exponential');

    error = sqrt((pd_true.mu-pd.mu)^2);
    tr = tr +error;
    error_list(i-19) = error/dist;
    %y = pdf(pd, x_pdf);
    y = pdf(pd, x_pdf);
    y2 = pdf(pd_true, x_pdf);
    ss = 0;
    for w=1:length(y)
        ss = ss + (y2(w)-y(w))^2;
    end
    rmse = sqrt(ss/length(y));
    rmse_list(i-19) = rmse;
end
figure
plot([20:1000],error_list)
title("Exponential Distribution: Loss across Samples")
xlabel('Number of Samples') 
ylabel('Percetn Loss') 
disp("average error:" + tr/981)
disp("percent Error:" + (tr/981)/dist)

figure
plot([20:1000],rmse_list)
title("Exponential Distribution: Loss across Samples")
xlabel('Number of Samples') 
ylabel('RMSE Loss') 

