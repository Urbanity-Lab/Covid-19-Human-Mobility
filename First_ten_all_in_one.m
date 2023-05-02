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
%   for i = column:-1:2
%       cases_actual(:,i-1) = cum_case(:,i)-cum_case(:,(i-1));
%   end
cases_actual = max(cum_case,0);
%csvwrite('Actual_Case.csv',movmean(cases_actual,7,2))
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
    u_c(k,:) = 200*sum/5;
    %u_c(k,:) = 20*sum/5;
    %u_c(k,:) = 200
    %u_c(k,:) = 200*u(i+3,:)-200*u(i+3,:)+200-200;
    k = k + 1;
end
u_c = u_c(1:end,:);

%% Hankel DMD
% %% Creating X,X',U and U' matrices
 observation_st = 21;
 observation_end = 85; 
 em = 29;
 predictionwindow = 21;
 X = MA(:,observation_st:observation_end); %Points for creating dynamics
 X = myhenkel(X,em);
 Xp  = MA(:,observation_st+1:observation_end+1);
 Xp = myhenkel(Xp,em);
 Ups = u_c(:,observation_st-14:observation_end-14);
 %Ups = u_c(:,observation_st:observation_end);
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
 %[rowA colA] = size(A)
%   for i = 1:rowA
%      for j = 1:colA
%          if abs(A(i,j)) < 0.00005
%             A(i,j)= 0;
%          else
%             A(i,j) = A(i,j);
%      end
%      end
%   end
 figure(4)
 subplot(2,1,1)
 imagesc(A)
 title('A-Matrix')
 colorbar
 subplot(2,1,2)
 imagesc(B)
 title('B-Matrix')
 colorbar 
 [W,eigs] = eig(A);
 figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
axis([-1.1 1.1 -1.1 1.1]);

%  figure(2)
%  imagesc(B)
%  colorbar
% [rowA colA] = size(A)
% A_new = zeros(rowA, colA); 

 %% 
 X_K = Xp(:,end);
 %Pred = zeros(m,predictionwindow);
 Update = MA(:,observation_st:observation_end);
 [m n] = size(MA);
 %% 
 %X_K = Xp(:,end); % First Input for Prediction
 for i= 1:predictionwindow
      x_k = A(end-(m-1):end,:)*X_K + B(end-(m-1):end,:)*U_1(:,observation_end+1+i-14);
      %x_k = A(end-(m-1):end,:)*X_K + B(end-(m-1):end,:)*U_1(:,observation_end+1+i)
      Pred(:,i) = x_k; % 
      Update = [Update(:,2:end) x_k]; % Updating 
      Update_1 = myhenkel(Update,em);
      X_K = Update_1(:,end);     
  end
 %%
 Actual = MA(:,observation_end+1:observation_end+predictionwindow);
 %% Plot 
 %csvwrite('Actual_Cases_in_window.csv',Actual)
 %csvwrite('Predicted_Cases_in_window.csv',Pred)
 j = 1; 
 %for i = [1 2 3 8 17 18]
 %for i = [1 3 8 18]
 for i = 1:20  
     filename = ["St.Lucie","Marion",'Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade'];
     figure(1)
     %subplot(3,2,j)
     %subplot(4,5,j)
     %subplot(2,2,j)
     subplot(5,4,j)
     plot(1:predictionwindow,Actual(i,:),'k-+','linewidth',2)
     title(sprintf('County: %s', filename(21-i)))
     %title(sprintf('County: %s', filename(i)))
     hold on
     %subplot(3,2,j)
     %subplot(2,2,j)
     %subplot(4,5,j)
     subplot(5,4,j)
     plot(1:predictionwindow,Pred(i,:),'r-*','linewidth',2)
     xticks([1 7 14 21 28])
     xticklabels({'06/25/20','07/02/20','07/08/20','07/14/20','07/21/20'})
     xlabel('Prediction Window')
     ylabel('Cumulative Cases')
     lg = legend('Actual','Predicted');
     lg.Location = 'Northwest';
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
