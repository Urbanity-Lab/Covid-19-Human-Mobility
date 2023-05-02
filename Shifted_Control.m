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
% %% Visulization
% figure
% imagesc(mobility_data)
% title('Google Mobility Report')
% colorbar
% figure
% imagesc(cum_case)
% title('Cumulative case')
% colorbar
% figure
% imagesc(cases_actual)
% title('Actual Cases')
% colorbar
 figure
 bar(flip(pop))
 xlabel('Counties')
 xticks([1:1:20])
 xticklabels({'St. Lucie','Marion','Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade'})
 ylabel('Population')
 title('County Population')
 figure
 imagesc(flip(cases_pk))
 title('Case Per Hundred Thousands')
 xlabel('Days')
 yticks([1:1:20])
 yticklabels({'St. Lucie','Marion','Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade'})
 title('Case Per Hundred Thousands')
 colorbar
 figure
 imagesc(flip(MA))
 xlabel('Days')
 yticks([1:1:20])
 yticklabels({'St. Lucie','Marion','Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade'})
 title('County Population')
 colorbar
 title('7 Day Moving Average of Case Per Hundred Thousands')
 colorbar

%% Control Data
maxi = max(mobility_data, [], 2); 
mini = min(mobility_data, [], 2);
lim = maxi-mini;
 for i = 1:column
     u(:,i) = (mobility_data(:,i)-mini)./lim;
 end
 k = 1;
 for i = 1:5:96
    sum = u(i,:);
    for j = 1:4
        sum = sum + u(i+j,:);
    end
    u_c(k,:) = 120*sum/5;
    k = k + 1;
end

figure
imagesc(u_c)
title('Control Input')
colorbar

%% Hankel DMD
% %% Creating X,X',U and U' matrices
 observation_st = 21;
 observation_end = 80; 
 em = 13;
 predictionwindow = 40;
 X = MA(:,observation_st:observation_end); %Points for creating dynamics
 X = myhenkel(X,em);
 Xp  = MA(:,observation_st+1:observation_end+1);
 Xp = myhenkel(Xp,em);
 Ups = u_c(:,observation_st-14:observation_end-14);
 Ups = myhenkel(Ups,em);
 U_1 = u_c(:,:);
 U_1 = myhenkel(U_1,em);
 %% SVD Analysis 
 Omega = [X;Ups];
 [U,Sig,V] = svd(Omega,'econ');
 thresh = 1e-9;
 rtil = length(find(diag(Sig)>thresh));
 U    = U(:,1:rtil); 
 Sig  = Sig(1:rtil,1:rtil);
 V    = V(:,1:rtil);
 %% A and B matrix calculation 
 A = Xp(1:end,:)*V*inv(Sig)*U(1:end/2,1:end)';
 B = Xp(1:end,:)*V*inv(Sig)*U(end/2+1:end,1:end)';
 %% 
 X_K = Xp(:,end);
 %Pred = zeros(m,predictionwindow);
 Update = MA(:,observation_st:observation_end);
 [m n] = size(MA);
 %% 
 %X_K = Xp(:,end); % First Input for Prediction
 for i= 1:predictionwindow
      x_k = A(end-(m-1):end,:)*X_K + B(end-(m-1):end,:)*U_1(:,observation_end+1+i-14);
      Pred(:,i) = x_k; % 
      Update = [Update(:,2:end) x_k]; % Updating 
      Update_1 = myhenkel(Update,em);
      X_K = Update_1(:,end);     
  end
 %%
 Actual = MA(:,observation_end+1:observation_end+predictionwindow);
 %% Plot 
 Diff = (Actual-Pred)./Actual*100;
 figure
 imagesc(flip(Actual));
 xlabel('days')
 ylabel('Counties')
 yticks([1:1:20])
 yticklabels({'St. Lucie','Marion','Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade'})
 title('Actual Case Count')
 colorbar
 figure
 imagesc(flip(Pred));
 xlabel('days')
 ylabel('Counties')
 yticks([1:1:20])
 yticklabels({'St. Lucie','Marion','Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade'})
 title('Prdicted')
 colorbar
 figure
 imagesc(flip(Diff(:,:)));
 xlabel('days')
 ylabel('Counties')
 yticks([1:1:20])
 yticklabels({'St. Lucie','Marion','Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade'})
 title('Difference')
 colorbar
 MSE  = mean(Diff(:,:).^2);
 RMSE = sqrt(MSE)
 figure
 plot(RMSE)
 %% Hankel Function.
 function hm = myhenkel(A,l)
     [m,n] = size(A);
     hm = zeros(m*(l+1), n-l);
     for k = 1:l+1
         hm(m*(k-1)+1:m*k,1:n-l) = A(1:m,k:k+n-l-1);
     end
 end
