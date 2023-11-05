% -----------------------------------------------------------------------
% JOINT OPTIMIZATION OF COMPUTATION OFFLOADING 
% AND RESOURCE ALLOCATION IN C-RAN WITH MOBILE EDGE
% COMPUTING USING EVOLUTIONARY ALGORITHM 
%------------------------------------------------------------------------
% Author - Sumit Singh
% NOTE- you MUST INSTALL CVX before running this code!
% This code is used to generate the dataset which is then used
% to train the a regression model. Model training has been done using
% MATLAB's APP.
%%
clear all; close all;
tic % measuring start time
Fontsize = 12;
global u genetic_fast e omega ohm gamma_C gamma_T Dm Fm B F N M Ln fm_local Tm_max Wm Sm kappa omega pf pt qb qc trainedModel3
% data_set ={};
iter  = 10; % ending UE count 10+iter*delta
start = 2; % starting UE count 10+start*delta 
delta = 10; % step size for incrementing number of UEs 
N = 10; %number of RRHs
load trainedModel3.mat 
%% creating random distribution of RRHs and UEs
cov_rad = 500;
% path loss model 37.6*log(dist) + 148.1
sigma_shadow = 8; %small scale fading model is independently and identically distributed (i.i.d.) Rayleigh fading with zero mean and unit variance
noise_power = -174;
xx0=0; yy0=0; %centre of disk
areaTotal = pi*cov_rad^2; %area of disk
% theta = 2*pi*(rand(N,1));
% rho = cov_rad*sqrt(rand(N,1));
theta = [asin(1/3)*ones(3,1);zeros(4,1);asin(-1/3)*ones(3,1)];
rho = cov_rad;
% [xxRRH, yyRRH] = pol2cart(theta,rho);
% xxRRH = xxRRH + xx0;
% yyRRH = yyRRH + yy0;
xxRRH = [-cov_rad/2,0,cov_rad/2,-3*cov_rad/4,-cov_rad/4,cov_rad/4,3*cov_rad/4,-cov_rad/2,0,cov_rad/2]' + xx0;
yyRRH = [1*cov_rad/2*ones(3,1);zeros(4,1);-1*cov_rad/2*ones(3,1)] + yy0;


 for q = start:iter
M = 10+delta*q; % number of UEs

thetaUE = 2*pi*(rand(M,1));
rhoUE = cov_rad*sqrt(rand(M,1));

[xxUE, yyUE] = pol2cart(thetaUE,rhoUE);
xxUE = xxUE + xx0;
yyUE = yyUE + yy0;

distance_matrix = pdist2([xxUE, yyUE],[xxRRH, yyRRH])/1000; % for distance in km
pathloss = 148.1+37.6*log10(distance_matrix) + 8*randn(1) - 10;
pt = 10.^30*ones(1,N); %30dbm
variance = 10.^(-pathloss/10); %gt
noise_variance = 10^-174;
interference = zeros(M,N);
for m = 1:M
    for n =1:N
        for j = 1:N
            if j~=n
                interference(m,n) = interference(m,n)+ pt(n)*variance(m,j);
            end
        end
    end
end
           
for m = 1:M
    for n =1:N
        e(m,n) = log2(1 + pt(n)*variance(m,n)/(noise_variance + interference(m,n)) );
    end
end
%%
W = zeros(M,3); %task matrix, each row cooresponds to one user's (Fm,Dm and Tm,max)
B = 1e7;    % channel bandwidth 100MHz
Ln = 5e7; % fronthaul capacity for each RRH is 50Mbps

F = 100e9; %MEC server maximum computational capacity 100 GHz
fm_local = 0.7e9;   % local computational capacity 0.7GHz
Dm = unifrnd(50,200,M,1)*1024*8; % output data size
Tm_max = unifrnd(0.6,1.5,M,1)*600e-3; % latency constraint of 600ms
Fm = 1000*Dm; % amountof computation for task Fm
Sm = 0.1; % storage cost 0.1$
kappa = .5;   % impact factor of computation cost
omega = 10;   % impact factor of badnwidth cost
pf = 0.03e-6; % unit price of charge for computation 0.03$/Mega Cycle
pt = 0.3/(1024*1024); %unit price of charge for transmission 0.3$
qb = 0.5e-6;   %unit price of charge for bandwidth 0.5$/Mhz
qc = 0.005e-6; %unit cost for computation 0.005$/mega cycle/s
Wm = [Fm,Dm,Tm_max]; % taks matrix
gamma_T = omega*B*qb*ones(M,1); % communication cost vector
gamma_C = kappa*F*qc*ones(M,1); % computationa cost vector
ohm = (pf*Fm +pt.*Dm+Sm); % revenue vector



[e_m,idx]= max(e,[],2);
a= zeros(M,N);

for m=1:M
    a(m,idx(m)) = 1;   
end

for i = 1:5
    permitted_ues = randi([0 1],M,1);
    [profit(i),sol(i),latency_check,~] = get_profit(permitted_ues.*a);
    data_set{q-1,i}=[nnz(permitted_ues),sol(i).N_new,profit(i),sol(i).b_sum,sum(sol(i).cm),latency_check,sol(i).ues_served,sol(i).mean_e,sol(i).total_outdata,sol(i).b_error_rmse,mean(mean(sol(i).b)),sol(i).cm_error_rmse,mean(sol(i).cm)];
end

 end
%  reshaped_set = reshape(data_set,prod(size(data_set)),1);
table = cell2table(vec(data_set));
table2 = splitvars(table);
table2.Properties.VariableNames = {'permitted_ues','N','profit','frac_RRH1','frac_RRH2','frac_RRH3','frac_RRH4','frac_RRH5','frac_RRH6','frac_RRH7','frac_RRH8',...
    'frac_RRH9','frac_RRH10','total_fac_cm','latency_check','num_RRH1','num_RRH2','num_RRH3','num_RRH4','num_RRH5','num_RRH6','num_RRH7','num_RRH8','num_RRH9','num_RRH10' ...
    ,'eff_RRH1','eff_RRH2','eff_RRH3','eff_RRH4','eff_RRH5','eff_RRH6','eff_RRH7','eff_RRH8','eff_RRH9','eff_RRH10','total_out_data','b_error_rmse','b_mean','cm_error_rmse','cm_mean'};
% table = cell2table(vec(data_set),'VariableNames',{'Permitted_UEs','Number of RRH','Profit','% of B @RRH 1',...
%    '% of B @RRH 2','% of B @RRH 3','% of B @RRH 4','% of B @RRH 5','% of B @RRH 6','% of B @RRH 7',...
%   '% of B @RRH 8','% of B @RRH 9','% of B @RRH 10','% of total Comput. Res.','Tmax-(T_Tr+Texe)'});

T = table2;
Q = T(T.profit>0 & T.latency_check >0,:); % for removing infeasible solutions
S = Q(:,logical([ones(1,14),0,[ones(1,21)]]));
% writetable(T,'raw_dataset.xlsx');
% writetable(Q,'filtered_dataset.xlsx');
writetable(T,'raw_dataset_error.xlsx');
writetable(Q,'filtered_dataset_error.xlsx');