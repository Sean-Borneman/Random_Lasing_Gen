
%Fit probability Distribution to data

load examgrades
x = grades(:, 1);
x_pdf = [1:0.1:100];
num_samples = 100;
x_test = reshape([1:1:num_samples], [], 1);

%---NORMAL------
disp("Exponential")
mu_true = 40;
sigma_true = 10;
pd_true = makedist('Exponential');
y_true = pdf(pd_true, x_test);
%rng('default') % For reproducibility
for i = 1:num_samples
    y_true(i) = random(pd_true);
end

pd = fitdist(y_true, 'Exponential');
disp("mu_true:" + pd_true.mu)

y = pdf(pd, x_pdf);
y2 = pdf(pd_true, x_pdf);
figure('Name', "Exponential fitting")
%histogram(y_true,'Normalization','pdf');
line(x_pdf, y2-y);
title("Normal Loss in Exponential Fitting")
xlabel('Range') 
ylabel('Probability') 

figure('Name', "Exponential fitting")
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
figure('Name', "Gamma fittin")
%histogram(y_true,'Normalization','pdf')
line(x_pdf, y2-y);
title("Loss in Gamma Fitting")
xlabel('Range') 
ylabel('Probability') 



figure('Name', "Gamma fitting")
%histogram(y_true,'Normalization','pdf');
line(x_pdf, y2);
line(x_pdf, y, 'Color', 'red','LineStyle','--');
legend({'Loss','Sampling Distribution'},'Location','northeast')
title("Gamma True V. Fitt PDF")
xlabel('Range') 
ylabel('Probability') 
