close all
clear all
clc
T = readtable('RawData110[uJcm2].csv');
F = readmatrix('Frequencies.csv');
T = table2array(T);
sample_num = size(T)
fake_num = 100
%x_pdf = F(1,:)
x_pdf = [100:1:500]
%x_pdf = x_pdf(1,:)
probs = zeros(1340,100);
for sample = 1:sample_num(1)
    pulse = T(:,sample );
    
    %---Kernel------
    pd = fitdist(pulse, 'Kernel','Bandwidth', 1);
    for i=1:fake_num
        probs(sample,i) = random(pd);
    end
    y = pdf(pd, x_pdf);
    figure('Name', 'Kernel fit')
    title("Kernel Distribution Fit")
    histogram(pulse,'Normalization','pdf')
    line(x_pdf, y)
    title("Fit Kernel Distribution to Data")
    xlabel('Range') 
    ylabel('Probability') 
end
imagesc(probs)