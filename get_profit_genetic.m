

function [profit_out,sol,x0]=get_profit_genetic(permitted_ues,a_se,e_global)
% A function to calculate profit for the genetic algorithm
global genetic_fast N Wm B  F Ln fm_local Sm kappa omega pf pt qb qc 
clear gamma_T

if genetic_fast ~= 1 % if normal genetic algorithm is used
idx = permitted_ues;
a = idx'.*a_se; % create offloading strtegy based on permitted UEs
profit_out = genetic_cvx(a); % get profit using cvx optimization solver
if isnan(profit_out)
    profit_out = 5*rand;
end
fprintf('profit_out = %f\n',profit_out);
else
%% Code for fast genetic algorithm
% Sm = 0.1; % storage cost 0.1$
% kappa = .5; % impact factor of abndwidth cost
% omega = 10; % impact factor of computation cost
% pf = 0.03e-6 ;% unit price of charge for computation 0.03$/Mega Cycle
% pt = 0.3/(1024*1024); %unit price of charge for transmission 0.3$
% qb = 0.5e-6; %unit price of charge for bandwidth 0.5$/Mhz
% qc = 0.005e-6;%unit cost for computation 0.005$/mega cycle/s



idx = permitted_ues;
disp(idx)
a = idx'.*a_se;
 
Wm_genetic = Wm(idx,:);
e = e_global(idx,:);
M = nnz(permitted_ues);
Dm = Wm_genetic(:,2);
Fm = Wm_genetic(:,1);
Tm_max =Wm_genetic(:,3);


cm = Fm./sum(Fm);   % computational resource allocation assumption
b = (0.2 +.3*rand)*(a)./sum(a); % radio resource allocation assumption
b(isnan(b)) = 0;

gamma_T = omega*B*qb*ones(M,1); %communication cost
gamma_C = kappa*F*qc*ones(M,1); % computation cost
ohm = (pf*Fm +pt*Dm+Sm); % revenue calculation
R = a.*b.*e*B; % transmission rates between UE and RRH
T_tr = Dm./sum(R')'; % transmission time 
T_exe = Fm./(cm*F); % execution time


% wrting the constraints
c1 = sum(Tm_max -(T_tr+T_exe)); % sum of all latency constraint violations
c2 = sum((fm_local*ones(M,1)-cm*F)*1e-10); % sum of violations of local computation resource less the cloud resource
c3 = sum((sum(R) - Ln*ones(1,N))*1e-8); % sum of violations of transmission rate available for a UE vs fronthaul link capapcity of RRH

fprintf('c1 =%f, c2 =%f, c3 =%f \n',c1,c2,c3);

profit_out = sum(ohm - cm.*gamma_C - sum(b')'.*gamma_T) -100*(min(c1,0))^2 -10*(max(c2,0))^2 -10*(max(c3,0))^2; % the objective function for profit. The penalty terms for constraints violations are added 
if profit_out < 0
    profit_out = 1/(abs(profit_out));
end
fprintf('profit_out = %f \n',profit_out);
% disp(sum(b))
% disp(sum(cm))
%%
end

sol =0;
x0=0;
end




function [profit] = genetic_cvx(a)
global    e ohm gamma_C gamma_T Dm Fm B F   Ln fm_local Tm_max 

idx2keep_columns = sum(abs(a),1)>0 ;
idx2keep_rows    = sum(abs(a),2)>0 ;
a_new = a(idx2keep_rows,idx2keep_columns);
e_new = e(idx2keep_rows,idx2keep_columns);
Dm_new = Dm(idx2keep_rows);
Fm_new = Fm(idx2keep_rows);
M = nnz(idx2keep_rows);
N = nnz(idx2keep_columns);
ohm_new = ohm(idx2keep_rows);
gamma_C_new = gamma_C(idx2keep_rows);
gamma_T_new = gamma_T(idx2keep_rows);
Tm_max_new = Tm_max(idx2keep_rows);


% change the global variable s to local not use now

      
    
% cvx_solver sedumi
cvx_begin quiet
variables cm(M) b(M,N)  
expressions R(M,N) T_tr(M,1) T_exe(M,1) 

R = a_new.*b.*e_new*B;
T_tr = Dm_new.*inv_pos(sum(R')');
T_exe = Fm_new.*inv_pos((cm*F));
maximise(sum(ohm_new - cm.*gamma_C_new - sum(b')'.*gamma_T_new))
subject to
sum(R) <= Ln*ones(1,N);
0<=cm<=1;
sum(cm) <= 1;
sum(b) <= 1;
b>=0;
cm*F >= fm_local*ones(M,1);
T_tr+T_exe <= Tm_max_new;
cvx_end
profit = cvx_optval;

if (cvx_status ~= "Infeasible")
    feasible = 1; 
else
    profit =1;
end


end
