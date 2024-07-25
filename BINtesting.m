close all 
clear all
clc
%Fit probability Distribution to data

load examgrades
x = grades(:, 1);
x_pdf = [1:0.1:100];
num_samples = 20;
x_test = reshape([1:1:num_samples], [], 1);

%---NORMAL------
disp("Binomial")
N = 30;
p = 0.25;
dist = sqrt(N^2 + p^2);
pd_true = makedist('binomial', 'N', N, 'p', p);
y_true = pdf(pd_true, x_test);
rng('default') % For reproducibility
tr = 0;
error_list = [20:1:1000];
rmse_list = [20:1:1000];
for sn = 20:1000
    for i = 1:sn
        y_true(i) = random(pd_true);
    end

    pd = fitdist(y_true+abs(awgn(y_true,1)*0.01), 'binomial','NTrials',30);

    error = sqrt((pd_true.p-pd.p)^2);
    tr = tr +error;
    
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
disp("average error:" + tr/100)
disp("percent Error:" + (tr/100)/dist*100)

figure
plot([20:1000],rmse_list)
title("Binomial Distribution: Loss across Samples")
xlabel('Number of Samples') 
ylabel('RMSE Loss') 