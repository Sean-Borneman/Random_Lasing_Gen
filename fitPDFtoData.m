close all 
clear all
clc
%Fit probability Distribution to data

load examgrades
x = grades(:, 1);
x_pdf = [1:0.1:100];

%---NORMAL------
pd = fitdist(x, 'Normal');
m = mean(pd);
y = pdf(pd, x_pdf);
figure('Name', 'Normal Distribution fit')
title("Normal DIstribution Fit")
%histogram(x,'Normalization','pdf')
line(x_pdf, y)
title("FIt Normal Distribution to Data")
xlabel('Range') 
ylabel('Probability') 

%---BINOMIAL------
pd = fitdist(x, 'nbin');
m = mean(pd);
y = pdf(pd, x_pdf);
figure('Name', 'Binomial Distribution Fit fit')
title("Binomial Distribution Fit")
histogram(x,'Normalization','pdf')
line(x_pdf, y)
title("Binomial Distribution to Data")
xlabel('Range') 
ylabel('Probability') 

%---Kernel------
pd = fitdist(x, 'Kernel','Bandwidth', 1)
m = mean(pd);
y = pdf(pd, x_pdf);
figure('Name', 'Kernel fit')
title("Kernel Distribution Fit")
histogram(x,'Normalization','pdf')
line(x_pdf, y)
title("Fit Kernel Distribution to Data")
xlabel('Range') 
ylabel('Probability') 
%---Gama------
pd = fitdist(x, 'Gamma')
m = mean(pd);
y = pdf(pd, x_pdf);
figure('Name', 'Gamma fit')
histogram(x,'Normalization','pdf')
line(x_pdf, y)
title("Fit Gamma Distribution to Data")
xlabel('Range') 
ylabel('Probability') 
