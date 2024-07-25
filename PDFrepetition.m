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
disp("Normal")
mu_true = 40;
b = 10;
dist = sqrt(mu_true^2 + b^2);
pd_true = makedist('Normal','mu',mu_true,'sigma',b);
y_true = pdf(pd_true, x_test);
rng('default') % For reproducibility
tr = 0;
error_list = [20:1:1000];
size(error_list)

for sn = 20:1000
    y_true = reshape([1:1:sn], [], 1);
    for i = 1:sn
        y_true(i) = random(pd_true);
    end

    pd = fitdist(y_true, 'Normal');

    error = sqrt((pd_true.mu-pd.mu)^2+(pd_true.sigma-pd.sigma)^2);
    tr = tr +error;
    error_list(i-19) = error;
    %y = pdf(pd, x_pdf);
end
figure
plot([20:1000],error_list)
title("Normal Distribution: Loss across Samples")
xlabel('Number of Samples') 
ylabel('Loss') 
disp("average error:" + tr/981)
disp("percent Error:" + (tr/981)/dist*100)

%%--- Gamma------%%

a= 40;
b = 10;
dist = sqrt(a^2 + b^2);
pd_true = makedist('Gamma','a',a,'b',b);
y_true = pdf(pd_true, x_test);
rng('default') % For reproducibility
tr = 0;
for sn = 10:200
    for i = 1:sn
        y_true(i) = random(pd_true);
    end

    pd = fitdist(y_true, 'Gamma');

    error = sqrt((pd_true.a-pd.a)^2+(pd_true.b-pd.b)^2);
    if(error < 200)
        tr = tr +error;
    end
    y = pdf(pd, x_pdf);
end
disp("Gamma")
disp("average error:" + tr/100)
disp("percent Error:" + (tr/100)/dist*100)
