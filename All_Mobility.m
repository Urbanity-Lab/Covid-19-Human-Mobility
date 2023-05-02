clear all;
close all;
clc;
%% Data Collection
cum_case = readtable('Florida_Cumulative_3_30_to_7_31.csv');
mobility_data = readtable('Florida_Mobility_3_30_to_7_31.csv');
population = readtable('FL_county_population.csv');
%% Create Matrix From Table
cum_case = cum_case{:,:};
cum_case = cum_case(:,:);
mobility_data = mobility_data{:,:};
mobility_data = mobility_data(:,:)';
FL_population = population(:,:);
pop = FL_population{:,:};
pop = flip(sort(pop));
pop = pop(1:20,:);
[row column] = size(cum_case);
 for i = column:-1:2
     cases_actual(:,i-1) = cum_case(:,i)-cum_case(:,(i-1));
 end
cases_actual = max(cases_actual,0);
csvwrite('Actual_Case.csv',movmean(cases_actual,7,2))
cases_pk = 100000*cases_actual./pop;
MA = movmean(cases_pk,7,2);

%% Control Data

%  for i = 1:column
%      u(:,i) = mobility_data(:,i);
%  end
%  k = 1;
%  for i = 1:5:96
%     sum = u(i,:);
%     for j = 1:4
%         sum = sum + u(i+j,:);
%     end
%     u_c(k,:) = sum;
%     [m n] = size(u_c);
%     figure(1)
%     plot(1:n, movmean(u_c(k,:),7,2),'-o','linewidth',1.25)
%     title('Average County Mobility Change')
%     xticks([3 33 64 94 124])
%     xticklabels({'04/01/20','05/01/20','06/01/20','07/31/20'})
%     xlabel('Date')
%     ylabel('Change With Respect to Baseline')
%     legend('St.Lucie','Marion','Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade')
%     hold on
%     k = k + 1; 
%  end
% 
% for ini = 1:5 
%   sum_f = zeros(100,124);
%     for row = ini:5:96+ini
%     sum_f(ini,:) = sum_f(ini,:) + mobility_data(row,:);
%     end
%    sum_f(ini,:) = sum_f(ini,:)/20;
%    [m n] = size(sum_f);
%    figure(2)
%    plot(1:n, movmean(sum_f(ini,:),7,2),'-o','linewidth',1.5)
%    title('Florida Mobility Change (County Aggregated)')
%    xticks([3 33 64 94 124])
%    xticklabels({'04/01/20','05/01/20','06/01/20','07/31/20'})
%    xlabel('Date')
%    ylabel('Change With Respect to Baseline')
%    lg = legend('Retail and recreation' ,'Grocery and pharmacy', 'Transit stations', 'Parks', 'Work stations');
%    lg.Location = 'southeast';
%    lg.FontSize = 6;
%    hold on
% end 

for k = 1:20
    figure(3)
    plot(MA(k,:),'linewidth',1.5)
    xticks([3 18 33 48 64 79 94 110 124])
    xticklabels({'Apr 01','Apr 15','May 01','May 15','Jun 01','Jun 15','Jul 01','jul 15','jul 31'})
    xlabel('Timeline','Fontsize',14,'Color','b','Fontweight','bold')
    ylabel('New Cases Per 100k (7-day Moving Average)','Fontsize',14,'Color','b','Fontweight','bold') 
    hold on
    k = k + 1; 
 end



figure(3)
% x1 = xline(21,':r','linewidth',2.5);
% x1.Label = {'Traning Window Starts';'Apr 19'}
% x1.LabelVerticalAlignment = 'top';
% x1.FontSize = 12
% x1.FontWeight = 'bold'
% 
% x2 = xline(85,':r','linewidth',2.5);
% x2.Label = {'Traning Ends';'Jun 23'}
% x2.LabelVerticalAlignment = 'top';
% x2.FontSize = 12
% x2.FontWeight = 'bold'
% 
% x3 = xline(86,':b','linewidth',2.5);
% x3.Label = {'Prediction Starts';'Jun 24'}
% x3.LabelVerticalAlignment = 'top';
% x3.FontSize = 12
% x3.FontWeight = 'bold'
% 
% x4 = xline(114,':b','linewidth',2.5);
% x4.Label = {'Prediction Window Ends';'Jul 20'}
% x4.LabelVerticalAlignment = 'top';
% x4.FontSize = 12
% x4.FontWeight = 'bold'

v1 = [21 0; 85 0; 85 120; 21 120];
f1 = [1 2 3 4];
patch('Faces',f1,'Vertices',v1,'FaceColor','red','FaceAlpha',.09)
hold on
v2 = [85 0; 113 0; 113 120; 85 120;];
f2 = [1 2 3 4];
patch('Faces',f2,'Vertices',v2,'FaceColor','blue','FaceAlpha',.09)
hold on
text(43,110,'Training Window','Color','red','FontSize',16,'Fontweight','bold')
text(88,110,'Prediction Window','Color','blue','FontSize',16,'Fontweight','bold')

legend('St.Lucie','Marion','Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade')
hold on
