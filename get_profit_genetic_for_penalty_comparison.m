function [profit_out]=get_profit_genetic(permitted_ues,a_se,e_global)
% A function to calculate profit for the genetic algorithm
global genetic_fast N  Wm B  F Ln fm_local Sm kappa omega pf pt qb qc gen best best_feasible best_unpenalized lambda penalty_no trainedModel
   
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



idx = permitted_ues;
a = a_se(idx,:); 
Wm_genetic = Wm(idx,:);
e = e_global(idx,:);
M = nnz(permitted_ues);
N_new = nnz(sum(a,1)>0);

Dm = Wm_genetic(:,2);
Fm = Wm_genetic(:,1);
Tm_max =Wm_genetic(:,3);

% cm = Fm./sum(Fm);   % computational resource allocation assumption
% b = (0.2 +.3*rand)*(a)./sum(a); % radio resource allocation assumption
% b(isnan(b)) = 0;

ues_served= sum(a)'; % for using same name as in the regression model
M_for_pred = M*ones(10,1);
N_for_pred = N_new*ones(10,1);
T = table(M_for_pred,N_for_pred,ues_served,'VariableNames',{'permitted_ues','N','num_RRH1'});
% trainedModel.predictFcn(T)
% class(trainedModel.predictFcn(T))
b_perRRH = trainedModel.predictFcn(T).*logical(ues_served);
b= b_perRRH'.*(a)./ues_served';
b(isnan(b)) = 0;

% load cmModel.mat
% T2 = table(M,N_new,sum(Dm),'VariableNames',{'permitted_ues','N','total_out_data'});
% total_comp_frac = cmModel.predictFcn(T2);

 cm = 0.8*Fm./sum(Fm);   % computational resource allocation assumption

gamma_T = omega*B*qb*ones(M,1); %communication cost
gamma_C = kappa*F*qc*ones(M,1); % computation cost
ohm = (pf*Fm +pt*Dm+Sm); % revenue calculation
R = a.*b.*e*B; % transmission rates between UE and RRH
T_tr = Dm./sum(R')'; % transmission time 
T_exe = Fm./(cm*F); % execution time

objective = sum(ohm - cm.*gamma_C - sum(b')'.*gamma_T);

% wrting the constraints from optimization equation (13)
c3 = (T_tr+T_exe)- Tm_max;  %  all latency constraint violations
c4 = (fm_local*ones(M,1)-cm*F); %  local computation resource less the cloud resource
c5 = (sum(R) - Ln*ones(1,N)); %  transmission rate available for a UEs vs fronthaul link capapcity of RRHs
% fprintf('c3 =%f, c4 =%f, c5 =%f \n',c3,c4,c5);
switch penalty_no
    case 1
% orginal penalty function 
penalty = -100*(max(sum(c3),0))^2 -10*(max(sum(c4*1e-10),0))^2 -10*(max(sum(c5*1e-8),0))^2;

    case 2
% static penalty function
penalty = -30*sum(c3>0) -10*sum(c4>0) -20*sum(c5>0);

% static penalty function_2
    case 3
penalty = -50*sum(max(c3,0)) - 30*sum(max(c4*1e-10 ,0)) - 40*sum(max(c5*1e-8, 0)) ; % kappa =1 


% Dynamic penalty functions
    case 4
penalty = gen*(-50*sum(max(c3,0)) - 30*sum(max(c4*1e-10 ,0)) - 40*sum(max(c5*1e-8, 0)));
% penalty = -50*sum(max(c3,0))^2 - 30*sum(max(c4*1e-10 ,0))^2 - 40*sum(max(c5*1e-8, 0))^2 ; % kappa =1 

% adaptive penalty function

    case 5
penalty = lambda*(-3*sum(max(c3,0)) - sum(max(c4*1e-10 ,0)) - 2*sum(max(c5*1e-8, 0))) ;
% penalty = -100*sum(max(c3,0)) - 100*sum(max(c4*1e-10 ,0)) -100*sum(max(c5*1e-8, 0)) ; % kappa =1 


% adaptive penalty function NFT based

    case 6
NFTo = 20;
lambda_k = 2;
NFTi = NFTo/(1+lambda_k*gen);

    
% for penalty6 the formula is different than that given in paper.
penalty = -(max(best) - max(best_unpenalized))*((3*sum(min(Tm_max -(T_tr+T_exe),0)) + sum(max((fm_local*ones(M,1)-cm*F)*1e-10 ,0)) + 2*sum(max((sum(R) - Ln*ones(1,N))*1e-8, 0))))/NFTi;

end
% profit_out(1) = objective + penalty1; % the objective function for profit.
profit_out = objective + penalty;
% fprintf("\n obj= %f penalty1 = %f penalty2 = %f  penalty3 = %f penalty4 = %f penalty5 = %f penalty6 = %f\n",objective,penalty1,penalty2,penalty3,penalty4,penalty5,penalty6)



if penalty_no == 5 && profit_out > best(gen)
    best(gen) = profit_out;
    best_feasible(gen) = penalty > -10; %sum(Tm_max -(T_tr+T_exe) >0 & (fm_local*ones(M,1)-cm*F) <0 & (sum(R) - Ln*ones(1,N)) <0);
    best_unpenalized(gen) = objective;
end


if profit_out < 0
    profit_out = 1/(abs(profit_out));
end

%  fprintf('profit_out = %f \n',profit_out);

%%
end


end




function [profit] = genetic_cvx(a)
global    e ohm gamma_C gamma_T Dm Fm B F   Ln fm_local Tm_max 

idx2keep_columns = sum(a,1)>0 ; 
idx2keep_rows    = sum(a,2)>0 ;
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
