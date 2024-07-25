close all
clear all
clc
a = zeros(1, 1601);
size(a)
for i = 1: 1
    mu = 540;                                %// Mean
    sigma = 10;                            %// Standard deviation
    mu2 = 550;
    sigma2 = 1;
    step_size =0.25/4;
    y0 = ones(1,10*sigma/step_size+1)*0.005;
    %// Plot curve
    x = (-5 * sigma:step_size:5 * sigma) + mu;  %// Plotting range
    y1 = exp(- 0.5 * ((x - mu) / sigma) .^ 2) / (sigma * sqrt(2 * pi))*0.1;
    y2 = exp(- 0.5 * ((x - mu2) / sigma2) .^ 2) / (sigma2 * sqrt(2 * pi))*0.1;

    size(y1);
    size(y0);
    y = y2*0.5+y1+y0;
    y4 = awgn(y, 1);
    y = y+y4*0.001;
    %plot(x, y)
    axis([510 570 0 inf])
    a(i, :) = y;
end
plot(a(1, :))
csvwrite("GenData.csv", a)