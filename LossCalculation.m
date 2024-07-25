
%Fit probability Distribution to data
close all
clear all
clc
load examgrades
x = grades(:, 1);
x_pdf = [1:0.1:100];
num_samples = 1000;
x_test = reshape([1:1:num_samples], [], 1);

%---NORMAL------
disp("Normal")
mu_true = 40;
sigma_true = 10;
pd_true = makedist('Normal', 'mu', mu_true, 'sigma', sigma_true);
y_true = pdf(pd_true, x_test);
%rng('default') % For reproducibility
for i = 1:num_samples
    y_true(i) = random(pd_true);
end

pd = fitdist(y_true, 'Normal');
disp("mu_true:" + pd_true.mu)

y = pdf(pd, x_pdf);
y2 = pdf(pd_true, x_pdf);
ss = 0;
for i=1:length(y)
    ss = ss + (y2(i)-y(i))^2;
end
rmse = sqrt(ss/length(y));
disp("rmse" + rmse);
figure('Name', "Normal fitting")
%histogram(y_true,'Normalization','pdf');
line(x_pdf, y2-y);
title("Loss in Normal Fitting")
xlabel('Range') 
ylabel('Loss') 

figure('Name', "Normal fitting")
%histogram(y_true,'Normalization','pdf');
line(x_pdf, y2);
line(x_pdf, y, 'Color', 'red','LineStyle','--');
legend({'Loss','Sampling Distribution'},'Location','northeast')
title("Normal True V. Fitt PDF")
xlabel('Range') 
ylabel('Probability') 
%--binomial---
disp("Binomial")
pd_true = makedist('binomial', 'N', 30, 'p', 0.25);
for i = 1:num_samples
    y_true(i) = random(pd_true);
end
pd = fitdist(y_true, 'binomial','NTrials',30);
disp("p_true:" + pd_true.p)
disp("p_exp:" + pd.p)
y = pdf(pd, x_pdf);
y2 = pdf(pd_true, x_pdf);
ss = 0;
for i=1:length(y)
    ss = ss + (y2(i)-y(i))^2;
end
rmse = sqrt(ss/length(y));
disp("rmse" + rmse);
figure('Name', "Binomial fitting")
%histogram(y_true,'Normalization','pdf')
line(x_pdf, y2-y)
title("Loss in Binomial FItting")
xlabel('Range') 
ylabel('Loss') 


figure('Name', "Binomial fitting")
%histogram(y_true,'Normalization','pdf');
line(x_pdf, y2);
line(x_pdf, y, 'Color', 'red','LineStyle','--');
legend({'Loss','Sampling Distribution'},'Location','northeast')
title("Binomial True V. Fitt PDF")
xlabel('Range') 
ylabel('Probability') 

%-- Gamma---
disp("Gamma");
pd_true = makedist('Gamma', 'a', 30, 'b', 1);
for i = 1:num_samples
    y_true(i) = random(pd_true);
end 
pd = fitdist(y_true, 'Gamma');
disp("a_true:"+pd_true.a);
disp("a_exp:"+pd.a);
disp("b_true:"+pd_true.b);
disp("b_exp:"+pd.b);
y = pdf(pd, x_pdf);
y2 = pdf(pd_true, x_pdf);
ss = 0;
for i=1:length(y)
    ss = ss + (y2(i)-y(i))^2;
end
rmse = sqrt(ss/length(y));
disp("rmse" + rmse);
figure('Name', "Gamma fittin")
%histogram(y_true,'Normalization','pdf')
line(x_pdf, y2-y);
title("Loss in Gamma Fitting")
xlabel('Range') 
ylabel('Loss') 



figure('Name', "Gamma fitting")
%histogram(y_true,'Normalization','pdf');
line(x_pdf, y2);
line(x_pdf, y, 'Color', 'red','LineStyle','--');
legend({'Loss','Sampling Distribution'},'Location','northeast')
title("Gamma True V. Fitt PDF")
xlabel('Range') 
ylabel('Probability') 


%---Exp------
disp("Exponential")
mu_true = 40;
sigma_true = 10;
pd_true = makedist('Exponential', 'mu', mu_true);
y_true = pdf(pd_true, x_test);
%rng('default') % For reproducibility
for i = 1:num_samples
    y_true(i) = random(pd_true);
end

pd = fitdist(y_true, 'Exponential');
disp("mu_true:" + pd_true.mu)

y = pdf(pd, x_pdf);
y2 = pdf(pd_true, x_pdf);
ss = 0;
for i=1:length(y)
    ss = ss + (y2(i)-y(i))^2;
end
rmse = sqrt(ss/length(y));
disp("rmse" + rmse);
figure('Name', "Exponential fitting")
%histogram(y_true,'Normalization','pdf');
line(x_pdf, y2-y);
title("Loss in Exponential Fitting")
xlabel('Range') 
ylabel('Loss') 

figure('Name', "Exponential fitting")
%histogram(y_true,'Normalization','pdf');
line(x_pdf, y2);
line(x_pdf, y, 'Color', 'red','LineStyle','--');
legend({'Loss','Sampling Distribution'},'Location','northeast')
title("Exponential True V. Fitt PDF")
xlabel('Range') 
ylabel('Probability') 

