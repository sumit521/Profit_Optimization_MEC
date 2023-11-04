function [profit_out,sol,constraint_check,a, T_tr, T_exe,exitflag ] = get_profit(a)
% A function to calculate proft by solving optimization problem using
% MATLAB's inbuilt solver

global e ohm gamma_C gamma_T Dm Fm B F   Ln fm_local Tm_max trainedModel3

% Resize the matrices and vectores
idx2keep_columns = sum(a,1)>0 ; 
idx2keep_rows    = sum(a,2)>0 ;
a_new = a(idx2keep_rows,idx2keep_columns);
e_new = e(idx2keep_rows,idx2keep_columns);
Dm_new = Dm(idx2keep_rows);
Fm_new = Fm(idx2keep_rows);
M = nnz(idx2keep_rows);
N = nnz(idx2keep_columns);
fprintf('M = %d, N= %d',M,N);

profit_max = optimproblem('ObjectiveSense','maximize');
% a = optimvar('a',[M,N],'Type','integer','LowerBound',0,'UpperBound',1);
b = optimvar('b',[M,N],'Type','continuous','LowerBound',0,'UpperBound',1);
cm = optimvar('cm',[M,1],'Type','continuous','LowerBound',0,'UpperBound',1);



R = a_new.*b.*e_new*B;
T_tr = Dm_new./sum(R')';
T_exe = Fm_new./(cm*F);

ohm_new = ohm(idx2keep_rows);
gamma_C_new = gamma_C(idx2keep_rows);
gamma_T_new = gamma_T(idx2keep_rows);
Tm_max_new = Tm_max(idx2keep_rows);
% Objective function
profit = sum(ohm_new - cm.*gamma_C_new - sum(b')'.*gamma_T_new);



%constraints
profit_max.Constraints.profit = profit >=0; %profit should be positive
profit_max.Constraints.computation = sum(cm) <=1;
profit_max.Constraints.bandwidth = sum(b) <=1;
profit_max.Constraints.fronthaul = sum(R) <= Ln*ones(1,N);
profit_max.Constraints.better_than_local = cm*F >= fm_local*ones(M,1);
profit_max.Constraints.Latency = T_tr+T_exe  <= Tm_max_new; 
profit_max.Constraints.profit = profit >=0; % profit should be positive
%options = optimoptions('fmincon','Algorithm','sqp','MaxIterations',1500,'MaxFunctionEvaluations',2000);%'ConstraintTolerance',1e-6,

% Initial conditions
x0.b = 2*rand(M,N)/(M*N);
x0.cm = 2*rand(M,1)/M;

profit_max.Objective = profit; 
%[profit_out,sol] = use_cvx(a_new,e_new, Dm_new, Fm_new, M, N, ohm_new, gamma_C_new, gamma_T_new, Tm_max_new);
solved =0;
counter =0;
loop_count = 0;
while(solved~=1)
    try 
        [sol, profit_out,exitflag]= solve(profit_max,x0);%,'options',options);
        if ~(exitflag == 1 || exitflag == 2 || exitflag == 3)
        profit_out = 0;
%         [profit_out,sol] = use_cvx(a);
        disp("infeasible problem");
        end
    
    catch ME
        fprintf("%s\n",ME.message);
        x0.b = 2*rand(M,N)/(M*N);
        x0.cm = 2*rand(M,1)/M;
        counter = counter +1;
        solved = 0;
    end
    
    if(counter == 0)
        solved =1;
        disp(profit_out)
    end
    loop_count = loop_count+1;
    if(loop_count == 10)
         %[profit_out,sol] = use_cvx(a_new,e_new, Dm_new, Fm_new, M, N, ohm_new, gamma_C_new, gamma_T_new, Tm_max_new);
        solved =1;
        disp("The offloading matrix is creating numerical problem. Try different one!");
         profit_out = 0;
        sol =0;
    end
    
        
end


R = a_new.*sol.b.*e_new*B;
T_tr = Dm_new./sum(R')';
T_exe = Fm_new./(sol.cm*F);

 constraint_check = mean(Tm_max_new -(T_tr+T_exe) ); 
 constraint_check % for displaying constraint on command window
 sum(sol.b)
 sum(sol.cm)
 sol.ues_served = sum(a);

 temp = zeros(1,length(idx2keep_columns));
 temp(idx2keep_columns)= sum(sol.b);
 sol.b_sum = temp;
 sol.N_new = N;
 sol.M_new = M;
 sol.total_outdata = sum(Dm_new);
 temp = zeros(1,length(idx2keep_columns));
 temp(idx2keep_columns)= mean(e_new);
 sol.mean_e = temp;
 
 
 %% for error comparison with prediction model
% ues_served= sum(a_new)'; % for using same name as in the regression model
% M_for_pred = M*ones(N,1);
% N_new = nnz(sum(a_new,1)>0);
% N_for_pred = N_new*ones(N,1);
% eff_RRH_for_pred = sol.mean_e(1)*ones(N,1);
% size(a_new)
% % disp(M_for_pred);
% % disp(N_for_pred);
% % disp(ues_served);
% %sprintf('%d,%d,%d \n',M_for_pred,N_for_pred,ues_served)
% T = table(M_for_pred,N_for_pred,ues_served,eff_RRH_for_pred,'VariableNames',{'permitted_ues','N','num_RRH1','eff_RRH1'});
% % trainedModel.predictFcn(T)
% % class(trainedModel.predictFcn(T))
% b_perRRH = trainedModel3.predictFcn(T).*logical(ues_served);
% b_pred = b_perRRH'.*(a_new)./ues_served';
% b_pred(isnan(b_pred)) = 0;
% sol.b_pred = b_pred;
% % load cmModel.mat
% % T2 = table(M,N_new,sum(Dm),'VariableNames',{'permitted_ues','N','total_out_data'});
% % total_comp_frac = cmModel.predictFcn(T2);
% 
% cm_pred = 0.8*Fm_new./sum(Fm_new);   % computational resource allocation assumption
% sol.cm_pred = cm_pred;
%  
% sol.b_error_rmse =  sqrt(mean((sol.b(:) - sol.b_pred(:)).^2));
% sol.cm_error_rmse = sqrt(mean((sol.cm(:) - sol.cm_pred(:)).^2));

end    




% the use_cvx function is optional, both fmincon and cvx give the same result

function [profit_out,sol] = use_cvx(a_new,e_new, Dm_new, Fm_new, M, N, ohm_new, gamma_C_new, gamma_T_new, Tm_max_new)
global      B F  Ln fm_local  
% change the globa variable s to local not use now

checker = 0;
counter =0;
while(checker ~= 1)       
    
% cvx_solver sedumi
cvx_begin 
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

if (cvx_status == "Solved")
    checker =1;
    profit_out = cvx_optval;
end
counter = counter+1;
if (counter == 3)
    profit_out = 0;
    checker =1;
end

end
sol = struct;
sol.b =b;
sol.cm =cm;
end




















%% working version of get_profit



% 
% function [profit_out,sol,a, T_tr, T_exe,exitflag ]=get_profit(a)
% global e ohm gamma_C gamma_T Dm Fm B F N  Ln fm_local Tm_max 
% % [profit_out,sol] = use_cvx(a);
% M = length(a(:,1));
% 
% profit_max = optimproblem('ObjectiveSense','maximize');
% % a = optimvar('a',[M,N],'Type','integer','LowerBound',0,'UpperBound',1);
% b = optimvar('b',[M,N],'Type','continuous','LowerBound',0,'UpperBound',1);
% cm = optimvar('cm',[M,1],'Type','continuous','LowerBound',0,'UpperBound',1);
% 
% 
% 
% R = a.*b.*e*B;
% T_tr = Dm./sum(R')';
% T_exe = Fm./(cm*F);
% 
% 
% profit = sum(ohm - cm.*gamma_C - sum(b')'.*gamma_T);
% 
% 
% %constraints
% profit_max.Constraints.profit = profit >=0; %profit should be positive
% profit_max.Constraints.computation = sum(cm) <=1;
% profit_max.Constraints.bandwidth = sum(b) <=1;
% profit_max.Constraints.fronthaul = sum(R) <= Ln*ones(1,N);
% profit_max.Constraints.better_than_local = cm*F >= fm_local*ones(M,1);
% profit_max.Constraints.Latency = T_tr+T_exe <= Tm_max;
% options = optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance',1e-6,'MaxIterations',1500,'MaxFunctionEvaluations',2000);
% x0.b = 2*rand(M,N)/(M*N);
% x0.cm = 2*rand(M,1)/M;
% feasibility_check = 0;
% counter =0;
% % while feasibility_check ~= 1
% %     %initial point
% % x0.b = 2*rand(M,N)/(M*N);
% % x0.cm = 2*rand(M,1)/M;
% % try 
% %         [sol, profit_out,exitflag] = solve(profit_max,x0,'options',options);
% %         if exitflag == 1 || exitflag == 2 || exitflag == 3
% %             feasibility_check = 1;
% %             disp("feasible problem");
% %         end
% %   catch ME
% %         fprintf("%s\n",ME.message);
% %         x0.b = 2*rand(M,N)/(M*N);
% %         x0.cm = 2*rand(M,1)/M;
% % end
% % 
% % 
% % 
% % 
% % 
% % counter = counter+1;
% % if counter == 3    
% %     feasibility_check = 1;
% %     disp("infeasible problem");      
% % end
% % end
% 
% 
% 
% 
% %objective function
% profit = sum(ohm - cm.*gamma_C - sum(b')'.*gamma_T);
% 
% profit_max.Constraints.profit = profit >=0; %profit should be positive
% profit_max.Objective = profit;
% 
% solved =0;
% counter =0;
% loop_count = 0;
% while(solved~=1)
%     try 
%         [sol, profit_out,exitflag]= solve(profit_max,x0,'options',options);
%         if ~(exitflag == 1 || exitflag == 2 || exitflag == 3)
%         profit_out = 0;
% %         [profit_out,sol] = use_cvx(a);
%         disp("infeasible problem");
%         end
%     
%     catch ME
%         fprintf("%s\n",ME.message);
%         x0.b = 2*rand(M,N)/(M*N);
%         x0.cm = 2*rand(M,1)/M;
%         counter = counter +1;
%         solved = 0;
%     end
%     
%     if(counter == 0)
%         solved =1;
%         disp(profit_out)
%     end
%     loop_count = loop_count+1;
%     if(loop_count ==10)
%          [profit_out,sol] = use_cvx(a);
%         solved =1;
%         disp("The offloading matrix is creating numerical problem. Try different one!");
%         profit_out = 0;
%         sol =0;
%     end
%     
%         
% end
% 
%     
% end    
% 
% 
% 
% 
% 
% 
% function [profit_out,sol] = use_cvx(a)
% global e ohm gamma_C gamma_T Dm Fm B F N M Ln fm_local Tm_max 
% checker = 0;
% counter =0;
% while(checker ~= 1)       
%     
% % cvx_solver sedumi
% cvx_begin 
% variables cm(M) b(M,N)  
% expressions R(M,N) T_tr(M,1) T_exe(M,1) 
% 
% R = a.*b.*e*B;
% T_tr = Dm.*inv_pos(sum(R')');
% T_exe = Fm.*inv_pos((cm*F));
% maximise(sum(ohm - cm.*gamma_C - sum(b')'.*gamma_T))
% subject to
% sum(R) <= Ln*ones(1,N);
% 0<=cm<=1;
% sum(cm) <= 1;
% sum(b) <= 1;
% b>=0;
% cm*F >= fm_local*ones(M,1);
% T_tr+T_exe <= Tm_max;
% cvx_end
% % msg = lastwarn;
% % if(contains(msg,"Matrix is close to singular or badly scaled. Results may be inaccurate."))
% %      disp("The offloading matrix is creating numerical problem. Try different one!");
% %      disp(msg);
% %       profit_out = 0;
% %       checker =1;
% %       msg ='';
% % end
% 
% if (cvx_status == "Solved")
%     checker =1;
%     profit_out = cvx_optval;
% end
% counter = counter+1;
% if (counter == 3)
%     profit_out = 0;
%     checker =1;
% end
% 
% end
% sol = struct;
% sol.b =b;
% sol.cm =cm;
% end
