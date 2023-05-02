clear all;
close all;
clc;
%% Data Collection
cum_case = readtable('Florida_Cumulative_3_30_to_7_31.csv');
%% Create Matrix From Table
cum_case = cum_case{:,:};
cases_actual = max(cum_case,0);
cases_pk = cases_actual;
MA = movmean(cases_pk,7,2);
MA = MA(1:end,:);

%% Hankel DMD
% %% Creating X,X',U and U' matrices
observation_st = 21;
observation_end = 85;
em = 55;
predictionwindow = 28;
X = MA(:,observation_st:observation_end); %Points for creating dynamics
X = myhenkel(X,em);
Xp  = MA(:,observation_st+1:observation_end+1);
Xp = myhenkel(Xp,em);
%% SVD Analysis 
Omega = X;
[U,Sig,V] = svd(Omega,'econ');
thresh = 1e-9;
rtil = length(find(diag(Sig)>thresh));
U    = U(:,1:rtil); 
Sig  = Sig(1:rtil,1:rtil);
V    = V(:,1:rtil);
%% A calculation 
A = Xp(1:end,:)*V*inv(Sig)*U(1:end,1:end)';
figure(1)
imagesc(A)
title('A-Matrix')
colorbar
[W,eigs] = eig(A);
figure(2)
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
axis([-1.1 1.1 -1.1 1.1]);
%% 
X_K = Xp(:,end);
%Pred = zeros(m,predictionwindow);
Update = MA(:,observation_st:observation_end);
[m n] = size(MA);
 %% 
 %X_K = Xp(:,end); % First Input for Prediction
 for i= 1:predictionwindow
      x_k = A(end-(m-1):end,:)*X_K
      Pred(:,i) = x_k; % 
      Update = [Update(:,2:end) x_k]; % Updating 
      Update_1 = myhenkel(Update,em);
      X_K = Update_1(:,end);     
  end
 %%
 Actual = MA(:,observation_end+1:observation_end+predictionwindow);
 %% Plot 
 j = 1; 
 for i = 1:20  
     filename = ["St.Lucie","Marion",'Lake','Osceola','Collier','Manatee','Sarasota','Seminole','Volusia','Pasco','Brevard','Polk','Lee','Duval','Pinellas','Orange','Hillsborough','Palm','Broward','Miami-Dade'];
     figure(1)
     
     subplot(5,4,j)
     plot(1:predictionwindow,Actual(i,:),'k-+','linewidth',2)
     title(sprintf('County: %s', filename(21-i)))
     %title(sprintf('County: %s', filename(i)))
     hold on
     subplot(5,4,j)
     plot(1:predictionwindow,Pred(i,:),'r-*','linewidth',2)
     xticks([1 7 14 21 28])
     xticklabels({'06/25/20','07/02/20','07/08/20','07/14/20','07/21/20'})
     %xlabel('Prediction Window')
     
     ax = gca;
     ax.YAxis.Exponent = 2
     if i==1
     a = annotation(gcf,'textarrow',...
    [0 0], [0.7 0.9],...
    'String','Cumulative Confirmed Cases', 'HeadStyle', 'none', 'LineStyle', 'none',...
    'FontSize',16, 'color','k','FontName', ...
    'cambria math','FontWeight','bold', 'TextRotation',90);
     %text([-2 2],[16 -16],'A Simple Plot','Color','red','FontSize',14)
     end
     if i>16
         xlabel('Date','Fontsize',14,'Fontweight','bold')
     end
     if i==20
     lg = legend('Actual Number of Cases','Prediction');
     %lg.Location = 'Southoutside';
     lg.FontSize = 11;
     lg.Orientation = 'horizontal';
     legend('boxoff')
     grid on
     end
     MeanSquaredError(j) = mean((abs((Actual(i,:) - Pred(i,:))/Actual(i,:))))*100
     j = j+1;
 end
 mean(MeanSquaredError) 
% Day_actual = [Actual(:,2:end)-Actual(:,1:end-1)];
% Day_pred = [Pred(:,2:end)-Pred(:,1:end-1)]; 
% daily_error = abs(Day_actual-Day_pred)/Day_actual*100
% figure(2)
% plot(daily_error)
% mean(daily_error)
 %% Hankel Function.
 function hm = myhenkel(A,l)
     [m,n] = size(A);
     hm = zeros(m*(l+1), n-l);
     for k = 1:l+1
         hm(m*(k-1)+1:m*k,1:n-l) = A(1:m,k:k+n-l-1);
     end
 end
