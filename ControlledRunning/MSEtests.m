line(mixxFit2,mixyfit2,'LineWidth',2,'Color',[1.00,0.41,0.16]);
line(realxFit2,realyfit2,'LineWidth',2,'Color',[1.00,0.41,0.16]);

d1 = realyfit2;
d1_len = numel(realyfit2);
d2 = mixyfit2;
d2_len = numel(mixyfit2);

d1 = interp1(1:(d2_len / d1_len):d2_len,d1,1:d2_len,'linear','extrap');
%line(mixxFit2,d1)

RMSE = sqrt(mean(((mixyfit2 - d1) .^ 2)))
% if 6 ==6
%     disp("true");
% end
% 
% slope = 0
%     abs(newPulseCollection(f-1, i) - newPulseCollection(f-5, i))+10;
%     x = 1:10;
%     y = newPulseCollection(f-10:f-1, i);
%     y = reshape(y,1,10);
%     fit = polyfit(x,y , 1);
%     y1 = polyval(fit, x);
%     figure
%     plot(x,y,'o')
%     hold on
%     plot(x,y1)
%     hold off
%     error = y1(1)-y1(length(y1));  