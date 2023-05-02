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
    u_c(k,:) = 300*sum/5;
    %u_c(k,:) = 20*sum/5;
    %u_c(k,:) = 200
    %u_c(k,:) = 200*u(i+3,:)-200*u(i+3,:)+200-200;
    k = k + 1;
end
% u_c = u_c(1:end,:);
u_c = u_c(1:end,:);

%% Hankel DMD
% %% Creating X,X',U and U' matrices
 observation_st = 21;
 observation_end = 85;
 el = 1;
 for em = 13:6:57
 %for em = 9:55
 %em = 21;
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
      x_k = A(end-(m-1):end,:)*X_K + B(end-(m-1):end,:)*U_1(:,observation_end+1+i-14);
      %x_k = A(end-(m-1):end,:)*X_K + B(end-(m-1):end,:)*U_1(:,observation_end+1+i)
      Pred(:,i) = x_k; % 
      Update = [Update(:,2:end) x_k]; % Updating 
      Update_1 = myhenkel(Update,em);
      X_K = Update_1(:,end);     
  end
 %%
 Actual = MA(:,observation_end+1:observation_end+predictionwindow);
 for i = 1:20  
     MeanSquaredError(i) = mean((abs((Actual(i,:) - Pred(i,:))/Actual(i,:))))*100;
 end
Average_error(el) = mean(MeanSquaredError);

 figure(1)
 theta = (0:1:100)*2*pi/100;
 subplot (2,4,el)
 plot(cos(theta),sin(theta),'r--') % plot unit circle
 hold on, grid on
 scatter(real(diag(eigs)),imag(diag(eigs)),'ob','filled')
 %axis([-1.1 1.1 -1.1 1.1]);
 xlim([-1.5 1.5])
 ylim([-1.5 1.5])
 xlabel('Real Axis')
 ylabel('Imaginary Axis')
 %set(gca,'DataAspectRatio',[10 1 1])
 %axis equal
 title("h = " + (em+1) + ", Error = " + round(Average_error(el),2) + "%")
 subtitle("W = 300")
  el = el+1;
 end
 %E = movmean(Average_error,4,2)
 figure(2)
 plot(Average_error)
 %% Hankel Function.
 function hm = myhenkel(A,l)
     [m,n] = size(A);
     hm = zeros(m*(l+1), n-l);
     for k = 1:l+1
         hm(m*(k-1)+1:m*k,1:n-l) = A(1:m,k:k+n-l-1);
     end
 end
