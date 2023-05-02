clear all;
close all;
clc;
%% Data Collection
cum_case = readtable('Florida_Cumulative_3_30_to_7_31.csv');
mobility_data = readtable('Mob.csv');
population = readtable('FL_county_population.csv');
%% Create Matrix From Table
cum_case = cum_case{:,:};
mobility_data = mobility_data{:,:};
mobility_data = mobility_data(:,:)';
FL_population = population(:,:);
pop = FL_population{:,:};
pop = flip(sort(pop));
pop = pop(1:20,:);
[row column] = size(mobility_data)
cases_actual = max(cum_case,0);
cases_pk = 100000*cases_actual./pop;
cases_pk = cases_actual;
MA = movmean(cases_pk,7,2);
MA = MA(1:end,:);

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
    u_c(k,:) = 600*sum/5;
    %u_c(k,:) = 20*sum/5;
    %u_c(k,:) = 200
    %u_c(k,:) = 200*u(i+3,:)-200*u(i+3,:)+200-200;
    k = k + 1;
end
% u_c = u_c(1:end,:);
u_c = u_c(1:end,:);
csvwrite('mobility.csv',u_c)
%% Hankel DMD
% %% Creating X,X',U and U' matrices
 observation_st = 21;
 observation_end = 85;
 el = 1;
 %for em = 15:8:56
 %for em = 9:55
 em = 55;
 predictionwindow = 28;
 X = MA(:,observation_st:observation_end); %Points for creating dynamics
 X = myhenkel(X,em);
 Xp  = MA(:,observation_st+1:observation_end+1);
 Xp = myhenkel(Xp,em);
 %Ups = u_c(:,observation_st-14:observation_end-14);
 Ups = u_c(:,observation_st:observation_end);
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
 [W,eigs] = eig(A);

 %% 
 X_K = Xp(:,end);
 %Pred = zeros(m,predictionwindow);
 Update = MA(:,observation_st:observation_end);
 [m n] = size(MA);
 
 for i= 1:predictionwindow
      %x_k = A(end-(m-1):end,:)*X_K + B(end-(m-1):end,:)*U_1(:,observation_end+1+i-14);
      x_k = A(end-(m-1):end,:)*X_K + B(end-(m-1):end,:)*U_1(:,observation_end+1+i);
      Pred(:,i) = x_k; % 
      Update = [Update(:,2:end) x_k]; % Updating 
      Update_1 = myhenkel(Update,em);
      X_K = Update_1(:,end);     
 end

 Actual = MA(:,observation_end+1:observation_end+predictionwindow);
 %% Plot 
 j = 1;
 for i = 1:20  
     filename = ["St.Lucie","Marion",'Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade'];
     figure(1)
     %subplot(3,2,j)
     subplot(4,5,j)
     %subplot(2,2,j)
     plot(1:predictionwindow,Actual(i,:),'*k','linewidth',1.5)
     title(sprintf('County: %s', filename(21-i)))
     hold on
     %subplot(3,2,j)
     %subplot(2,2,j)
     subplot(4,5,j)
     plot(1:predictionwindow,Pred(i,:),'*r','linewidth',1.5)
     xticks([1 7 14 21 28])
     xticklabels({'06/25/20','07/02/20','07/08/20','07/14/20','07/21/20'})
     xlabel('Prediction Window')
     ylabel('Cumulative Cases')
     lg = legend('Actual','Predicted');
     lg.Location = 'southeast';
     lg.FontSize = 6;
     grid on
     MeanSquaredError(j) = mean((abs((Actual(i,:) - Pred(i,:))/Actual(i,:))))*100
     j = j+1;
 end
 %% Hankel Function.
 function hm = myhenkel(A,l)
     [m,n] = size(A);
     hm = zeros(m*(l+1), n-l);
     for k = 1:l+1
         hm(m*(k-1)+1:m*k,1:n-l) = A(1:m,k:k+n-l-1);
     end
 end
