

function [profit_out,sol,x0]=get_profit_sjoora(a,Wm_sorted,e_m_sorted,test_count)
global B N F Ln fm_local Sm kappa omega pf pt qb qc

% Almost same as get_profit function. The only difference is that the constraint of 'fraction of computational resource to be allocated should be less than 1' has been removed. 




e = e_m_sorted;
M = test_count;
Dm = Wm_sorted(:,2);
Fm = Wm_sorted(:,1);
Tm_max =Wm_sorted(:,3);



profit_max = optimproblem('ObjectiveSense','maximize');

b = optimvar('b',[M,N],'Type','continuous','LowerBound',0,'UpperBound',1);
cm = optimvar('cm',[M,1],'Type','continuous','LowerBound',0,'UpperBound',1);
options = optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance',1e-6,'MaxIterations',1500,'MaxFunctionEvaluations',2000);
gamma_T = omega*B*qb*ones(M,1);
gamma_C = kappa*F*qc*ones(M,1);
ohm = (pf*Fm +pt*Dm+Sm);


% [profit_out,sol] = use_cvx(a);

R = a.*b.*e*B;
T_tr = Dm./sum(R')';
T_exe = Fm./(cm*F);

profit = sum(ohm - cm.*gamma_C - sum(b')'.*gamma_T);
profit_max.Objective = profit;

profit_max.Constraints.profit = profit >=0;
% profit_max.Constraints.computation = sum(cm) <=1;
profit_max.Constraints.bandwidth = sum(b) <=1;
profit_max.Constraints.fronthaul = sum(R) <= Ln*ones(1,N);
profit_max.Constraints.better_than_local = cm*F >= fm_local*ones(M,1);
profit_max.Constraints.Latency = T_tr+T_exe == Tm_max;

x0.b = rand(M,N);
x0.cm = rand(M,1);

solved =0;
counter =0;
loop_count = 0;
while(solved~=1)
    try 
        [sol, profit_out]= solve(profit_max,x0,'options',options);
    catch ME
        fprintf("%s\n",ME.message);
        x0.b = rand(M,N);
        x0.cm = rand(M,1);
        counter = counter +1;
        solved = 0;
    end
    if(counter == 0)
        solved =1;
        disp(profit_out)
    end
    loop_count = loop_count+1;
    if(loop_count ==10)
%         [profit_out,sol] = use_cvx(a);
        solved =1;
        disp("The offloading matrix is creating numerical problem. Try different one!");
        profit_out = 0;
        sol =0;
    end
    
        
end
 


end

function [profit_out,sol] = use_cvx(a)
global e ohm gamma_C gamma_T Dm Fm B F N M Ln fm_local Tm_max 
checker = 0;
counter =0;
while(checker ~= 1)       
    
% cvx_solver sedumi
cvx_begin 
variables cm(M) b(M,N)  
expressions R(M,N) T_tr(M,1) T_exe(M,1) 

R = a.*b.*e*B;
T_tr = Dm.*inv_pos(sum(R')');
T_exe = Fm.*inv_pos((cm*F));
maximise(sum(ohm - cm.*gamma_C - sum(b')'.*gamma_T))
subject to
sum(R) <= Ln*ones(1,N);
0<=cm<=1;
sum(cm) <= 1;
sum(sum(b)) <= 1;
b>=0;
cm*F >= fm_local*ones(M,1);
T_tr+T_exe <= Tm_max;
cvx_end
msg = lastwarn;
if(contains(msg,"Matrix is close to singular or badly scaled. Results may be inaccurate."))
     disp("The offloading matrix is creating numerical problem. Try different one!");
     disp(msg);
      profit_out = 0;
      checker =1;
      msg ='';
end

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
