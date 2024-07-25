close all 
clear all
clc
%Fit probability Distribution to data

load examgrades
x = grades(:, 1);
x_pdf = [1:0.1:100];
num_samples = 1000;
x_test = reshape([1:1:num_samples], [], 1);

%---NORMAL------
disp("Normal")
mu_true = 40;
sigma_true = 10;
pd_true = makedist('Normal','mu',mu_true,'sigma',sigma_true);
y_true = pdf(pd_true, x_test);
%rng('default') % For reproducibility
for i = 1:num_samples
    y_true(i) = random(pd_true);
end

pd = fitdist(y_true, 'Normal');
disp("mu_true:" + pd_true.mu)
disp("mu_exp:" + pd.mu)
disp("sigma_true:" + pd_true.sigma)
disp("sigma_exp:"  +pd.sigma)
y = pdf(pd, x_pdf);
yt = pdf(pd_true, x_pdf);
figure('Name', "Normal fitting")
%histogram(y_true,'Normalization','pdf');
line(x_pdf, yt);
title("Normal Distribution")
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
figure('Name', "Binomial fitting")
histogram(y_true,'Normalization','pdf')
line(x_pdf, y)

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
figure('Name', "Gamma fittin")
histogram(y_true,'Normalization','pdf')
line(x_pdf, y)
