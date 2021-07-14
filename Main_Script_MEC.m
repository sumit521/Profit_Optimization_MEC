% -----------------------------------------------------------------------
% PROFIT OPTIMIZATION FOR MOBILE EDGE COMPUTING USING GENETIC ALGORITHM 
%------------------------------------------------------------------------
% Author - Sumit Singh
% NOTE- you MUST INSTALL CVX before running this code!
% Using this code, we try to evaluate the performance of three algorithms
% in terms of their profitability in deciding computation offloading
% strategy in mobile edge computing scenario. These three algorithms are
% 1. Normal Genetic 2. Fast Genetic 3. SJOORA scheme.
% For system model, simulation paramters and SJOORA algorithm, refer to https://ieeexplore.ieee.org/document/8736324
%
% The program execution can take upto 45 min for u=3 (3 iterations) on i7-1170 processor,
% RAM 16GB. More time for higher u. Don't set u=1:1. As the mean() function
% will avarage the single row. Keep u>1, for averaging along columns.
% Have patience!
% Copyright 2021, Sumit Singh, All rights reserved
%%
clear all; close all;
tic % measuring start time
Fontsize = 12;
global u genetic_fast e omega ohm gamma_C gamma_T Dm Fm B F N M Ln fm_local Tm_max Wm Sm kappa omega pf pt qb qc
iter  = 10; % ending UE count 10+iter*delta
start = 2; % starting UE count 10+start*delta 
delta = 10; % step size for incrementing number of UEs 
for u = 1:3 % number of iterations to average the result
%%
for q = start:iter
%%
N = 10; %number of RRHs
M = 10+delta*q; % number of UEs
W = zeros(M,3); %task matrix, each row cooresponds to one user's (Fm,Dm and Tm,max)
B = 1e7;    % channel bandwidth 100MHz
%b = rand(M,N); % percentage of radio spectrum allocated to MT m by RRH n
p = 10*rand(M,N); % transmit power from RRH n to MT m ; for each RRH Tx power is 30dBm
g = 50*rand(M,N); % channel gain from RRH n to MT m
sigma_square = 5; %AWGN  noise with zero mean -174dBm/Hz
Ln = 5e7; % fronthaul capacity for each RRH is 50Mbps

F = 100e9; %MEC server maximum computational capacity 100 GHz
fm_local = 0.7e9;   % local computational capacity 0.7GHz
Dm = unifrnd(50,200,M,1)*1024*8; % output data size
Tm_max = unifrnd(0.6,1.5,M,1)*600e-3; % latency constraint of 600ms
Fm = 1000*Dm; % amountof computation for task Fm
Sm = 0.1; % storage cost 0.1$
kappa = .5;   % impact factor of abndwidth cost
omega = 10;   % impact factor of computation cost
pf = 0.03e-6; % unit price of charge for computation 0.03$/Mega Cycle
pt = 0.3/(1024*1024); %unit price of charge for transmission 0.3$
qb = 0.5e-6;   %unit price of charge for bandwidth 0.5$/Mhz
qc = 0.005e-6; %unit cost for computation 0.005$/mega cycle/s
Wm = [Fm,Dm,Tm_max]; % taks matrix
gamma_T = omega*B*qb*ones(M,1); % communication cost vector
gamma_C = kappa*F*qc*ones(M,1); % computationa cost vector
ohm = (pf*Fm +pt.*Dm+Sm); % revenue vector

% creating spectrum efficiency matrix  
e = 5*rand(M,N).*rand(M,N);% spectrum efficiency between each UE nad RRH is in range 0 to 5.



%% Calculating the profits 

% genetic_fast = 1 ;
% [profit_genetic(u,q),~] = genetic_algo();

for s = 1:2 
    genetic_fast = s-1; %if genetic_fast == 1, then fast genetic algo is used otherwise normal genetic
[profit_genetic(u,q,s),~] = genetic_algo();
end

[profit_sjoora(u,q),S1] = sjoora(); % calculate profit using sjoora algorithm

%%
end
end
toc % measuring end time

%% Plotting Results 
x = 10+delta*(start:iter); % X axis of the curves, representing number of UEs 

%% for resource dependency trend
% markers = {'o','+','*','s','d','v','>','h'};
% % List a bunch of colors; like the markers, they 
% % will be selected circularly. 
% colors = {'b','c','r','g','k','m'};
% 
% % Same with line styles
% linestyle = {'-','--','-.',':'};
% % this function will do the circular selection
% % Example:  getprop(colors, 7) = 'b'
% getFirst = @(v)v{1}; 
% getprop = @(options, idx)getFirst(circshift(options,-idx+1));
% 
% figure;
% hold on
% for i = 1:2:u
%     y = profit_genetic(i,:);    
%     plot(x,y(start:iter),'Marker',getprop(markers,i),...
%         'color',getprop(colors,i));
%     text{i} = sprintf('Ln = %d Mbps',(20*i));     
% end
% legend(text{1:2:u})
% xlabel("Number of UEs",'FontSize',Fontsize)
% ylabel("Profit (in $)",'FontSize',Fontsize)
% title("Impact of Fronthaul Link capacity per RRH",'FontSize',Fontsize+1)
% grid on;
% hold off


%% for all three algorithms comparison



% Comparison of two genetic algos
% profit_genetic_avg = mean(profit_genetic(:,:,1));
% profit_fast_avg = mean(profit_genetic(:,:,2));
% hold on
% plot(x,profit_genetic_avg(start:iter),'-m*')
% plot(x,profit_fast_avg(start:iter),'-kd')
% hold off
% legend('Normal Genetic','Fast Genetic')
% xlabel("Number of UEs",'FontSize',Fontsize)
% ylabel("Profit (in $)",'FontSize',Fontsize)
% title("Comparison of Genetic Algorithms",'FontSize',Fontsize+1)
% grid on;

% Comparison of all three algorithms (including sjoora)
profit_genetic_avg = mean(profit_genetic(:,:,1));
profit_fast_avg = mean(profit_genetic(:,:,2));
profit_sjoora_avg = mean(profit_sjoora);
figure;
hold on
plot(x,profit_genetic_avg(start:iter),'-m*')
plot(x,profit_fast_avg(start:iter),'-kd')
plot(x,profit_sjoora_avg(start:iter),'-co')

xlabel("Number of UEs",'FontSize',Fontsize)
ylabel("Profit (in $)",'FontSize',Fontsize)
title("Comparison of offloading strategies",'FontSize',Fontsize+1)
grid on;

hold off
legend('Normal Genetic','Fast Genetic','Modified SJOORA')




%% functions defined for use in script above


%% SJOORA Algorithm - Refer Algorithm 2 of https://ieeexplore.ieee.org/document/8736324

 function [profit_sjoora,S1] = sjoora()
 disp('Solving using SJOORA')
global e M N Wm
S1 = []; %set of allowed MTs (mobile terminals)/UEs(User Equipment)
S2 = []; % set of disallowd UEs
S3 = 1:M; % test set for deciding to put UEs in S1 or S2
x = zeros(M,N); % seed for creating inital offloading matrix
K = numel(S3); % number of UEs in test set

for m=1:M
    e_m(m) = max(e(m,:)); % find the highest spectrum efficiency available to a UE from all the RRHs
    n = find(e(m,:)==e_m(m)); % find the corresponding RRH
    x(m,n) = 1;  % set 1 for that RRH
 
end
e_m=repmat(e_m',1,N); % repeat the columns of e_m to create the matrix mentioned in Algorithm in paper



cm_total =0; % total fraction of MEC server computational resource allocated
while(K~=0)
fprintf('K = %d \n',K); % display value of K
[~,idx_S3] =sort(e_m(S3,1),'Descend'); % sort the MTs in terms of spectrum effiency (SE) 
% sort the other papameter also so that they are in same order as MTs
idx = S3(idx_S3); % the  MTs in S3 sorted on the basis of their SE
x_sorted = x(idx,:); % sort the offloading matrix
e_m_sorted = e_m(idx,:);% sort the spectrum efficiency matrix
Wm_sorted = Wm(idx,:); % sort the tasks of UE 

test_count = ceil(K/2); % select half the number of MTs in S3 for testing the feasibility
a = x_sorted(1:test_count,1:N);% keep only the corresponding half component of the sorted matrix 

[~,sol,~] =get_profit_sjoora(a,Wm_sorted(1:test_count,:),e_m_sorted(1:test_count,:),test_count);% solve the optimization problem for the selected MTs
b_sorted = sol.b; % the optimal radio resource allocation ratios *b
cm_sorted = sol.cm; % the optimal computational resource allocation ratios *c
Um = total_mt_profit(Wm_sorted(1:test_count,:),b_sorted,cm_sorted); % calculate the profit for the test count MTs with optimal radio and computational resource allocation arrived in previous step


idx_temp = idx(1:test_count); % list of test MTs of test count

cm_total = sum(cm_sorted)+cm_total; % add the computational resource allocation factor *c of optimzation problem solved above to the total  

fprintf('sum(cm_sorted)= %f \n',sum(cm_sorted)); %for testing purpose

if(cm_total<=1)  % if the there are still resources to be allocated, i.e.  there is no overshoot of computational resources then
    S1 = [S1,idx_temp]; % add the MTs of test count (K/2) to the existing set of permitted MTs S1
    S3 = setdiff(S3,S1); % remove the MTs present in set S1 from set S3
else
    cm_total = cm_total - sum(cm_sorted); % the computational resources consumption corresponding to selected UEs only. Reversing the cm added this iteration of while loop as no UE put in S1
    
    index_min = find(Um ==min(Um)); % find the index for MT with minimum profit
    rrh_min = find(a(index_min,:)~=0); % find the RRH through which the above MT caan offload
    disallowed_mts = []; % create a temporary empty set of allowed MTs
    allowed_mts =[];
    for i = 1:test_count % for all the MTs in test test count
        if(find(a(i,:)~=0)==rrh_min && e_m_sorted(i,1)<=e_m_sorted(index_min,1)) % find the the MTs which share the same RRH as the one with min. profit MT and have spectrum efficincy less than or equal to          
                disallowed_mts = [disallowed_mts,i]; % add them to set of disallowed MTs
                
                     
        else
            allowed_mts =[allowed_mts,i];% otherwise they are allowed
            
        end
          
    end
    fprintf('disallowed_UEs index ');
    disp(disallowed_mts)
    
    S2 = [S2,idx_temp(disallowed_mts)]; % add the disallowed MTs to set S2
    S3 = setdiff(S3,S2); % remove the MTs of set S2 from set S3

    
end
  
   K = numel(S3); % revise the value of K as the number of MTs in test set S3
   
end

disp('selcted UEs')
disp('S1')
disp(S1)
fprintf('cm_total = %f \n',cm_total)

    a =x; % recreate the inital offloading matrix 
    feasible = 0; % marker for checking the feasibility of set S1 
    counter = 1; % set the counter to 1
    count_limit = length(S1); % check till all MTs in S1 are exhausted
    end_limit = 2; % descrement size
    a_best = ismember(1:M,S1);% a vector of size M having ones at indexes corresponding to allowed UEs, 0 otherwise
    
    while feasible ~= 1 && counter ~= count_limit      
    non_zero_idx = find(a_best~=0); % find the indexes for allowed MTs   
    selected_idx = randsample(non_zero_idx,length(non_zero_idx)-end_limit);% randomly select 2(step size) less MTs
    a_best(~ismember(1:length(a_best),selected_idx)) = 0; % set the two selected MTs indexes to 0 or drop them
    [feasible] = feasibility_check(a_best'.*a); % check the feasibility of the offloading strategy arrived by above steps
    
    if ~mod(counter,3) % if the counter has ran three consecutive iterations
        end_limit = end_limit+1;% then increase the dropping step size by 1
    end
    counter = counter +1;
 
    end
    
    % again we calculate profit using fmincon and check for the violation of any resource constraint 
    counter = 1;
    count_limit = length(a_best);
    end_limit = 2;
    [profit_sjoora,~,constraint_check] = get_profit(a_best'.*a); % this time get profit using get_profit() function
    % again keep reducing UEs from a_best while checking for constraint
    while (constraint_check < -1e-2 || profit_sjoora == 0) &&  counter ~= count_limit
        
        non_zero_idx = find(a_best~=0);    
    selected_idx = randsample(non_zero_idx,length(non_zero_idx)-end_limit);
    a_best(~ismember(1:length(a_best),selected_idx)) = 0;
    [profit_sjoora,~,constraint_check] = get_profit(a_best'.*a);
    if ~mod(counter,3)
        end_limit = end_limit+1;
    end
    counter = counter +1;
    
    end
    
    total_mts = 1:M;
    S1 = total_mts(a_best>0); % ouput the S1 as chosen MTs for offloading 


end




function [profit_out]= total_mt_profit(Wm_sorted,b,cm)
% function to calculate the profit once radio and computational resource
% distribution known. It has been used in SJOORA function
global B  F   Sm kappa omega pf pt qb qc
M = length(cm);
Dm = Wm_sorted(:,2);
Fm = Wm_sorted(:,1);
Tm_max =Wm_sorted(:,3);


gamma_T = omega*B*qb*ones(M,1);
gamma_C = kappa*F*qc*ones(M,1);
ohm = (pf*Fm +pt.*Dm+Sm);
profit_out = ohm - cm.*gamma_C - sum(b')'.*gamma_T;
end

%%






%----------------------------------------------------------
% GENETIC ALGORITHM
%----------------------------------------------------------
function [profit_genetic, a_best] = genetic_algo()
% function to select the best UEs for offloading
% Refer Algorithm1 in https://ieeexplore.ieee.org/document/6253581 

global  genetic_fast e N M Wm  

if genetic_fast == 1
    nind = 2000; % if fast genetic, set the number of individuals in the population high
    disp('Solving using Fast Genetic Algorithm')
else
    nind =15; % if normal genetic, set the the number of individuals in the population low  
    disp('Solving using Normal Genetic Algorithm')
    
end
ggap = 0.8; % generation gap
mutr = 0.05; % mutation ratio
maxgen = 3; % maximum number of generations 
a = zeros(M,N);
% create inital offloading strategy matrix 
for m=1:M
    e_m(m) = max(e(m,:));
    n = find(e(m,:)==e_m(m));
    a(m,n) = 1;   
end
Wm_metric = [Wm(:,1)*1e-9,Wm(:,2)*1e-6,Wm(:,3)*1.5]; % create task metric
x1 = max(sum(Wm_metric')') - min(sum(Wm_metric')'); %calculating range of Wm
x2 = max(e_m) - min(e_m); % calculate range of e_m

rating = (x1*e_m + x2*sum(Wm_metric'))/(x1+x2); % deciding rating based on the the ratio of ranges of Wm and e_m


[~, sorted_ues] = sort(rating,'descend'); % sort the UEs in descending order of their ratings

step_size = ceil(0.1*M); % select the stepsize
nones = 20;
check(1)=0;% first index of check vector set to 0
i=2;
feasible = -5; % set intial condition as infeasible
while feasible ~= 0 && nones <M
    selected_ues = ismember(1:M,sorted_ues(1:nones));
    check(i) = feasibility_check(selected_ues'.*a);
    if check(i) == 0 && check(i-1) == 0 % if two consecutive attempts are infeasible
        feasible = 0; % then set the value of feasible variable to 0
    end
   nones = nones+step_size; % increment the number of ones by step_size
end
nones = min(nones-step_size,M); % for removing the increment caused by last infeasible case


for j= 1:nind
permitted_ues = randsample(1:M,nones,true,rating); % select permitted UEs based on nones and probability of a UE being selected is proportional to its rating
chrom(:,:,j) = ismember(1:M,permitted_ues); % create a chromosome of length M  
end

gen =0;
while(gen < maxgen) % while the maximum number of generations not reached
  
    for i = 1:nind
      objv(i) = get_profit_genetic(chrom(:,:,i),a,e); % calculate the profit using the function. Create objctive vector of profits 
    end       
        
        fitnv = objv./norm(objv); % create normalized fitness vector
        
        selch =rouelttewheelselect(chrom,fitnv,ggap); % roulette wheel selection
        new_nind = floor(ggap*nind); % new individuals to be created for next generation
        i=1;
        while(i+1<new_nind )
            [selch(:,:,i),selch(:,:,i+2)] = crossover(selch(:,:,i),selch(:,:,i+1)); % crossover of selected chromosomes
            i = i+2;
        end
        
        for i=1:new_nind
            selch(:,:,i) = mutation(selch(:,:,i),mutr); % mutation of selected chromosomes s per mutation ratio
        end
        for i=1:new_nind
            objvselch(i) = get_profit_genetic(selch(:,:,i),a,e); % again created an objective vector of profits for selected chromosomes          
        end
        chrom = reinsert(chrom,selch,objv,objvselch); % reinsert the new individuals in the population by slecting best among old and new such that toatl count remains nind   
    
        gen = gen+1;  
    
end

    for i = 1:nind
        objv(i) = get_profit_genetic(chrom(:,:,i),a,e);% calculate the profits for the final population after 3 generations
    end
    
    [objv_sorted,idx] = sort(objv,'descend'); % sort the profits in descending order
    
if genetic_fast ~= 1    
    profit_genetic = objv_sorted(1); % if it is normal genetic algo, set the output profit as top of sorted objv 
    a_best = idx(1); % best chromosome corresponding to chromosome with highest profit
    else %condition for fast genetic algorithm
%%   fast genetic algorithm
for i = 1:5 % for choosing best of 5 profits
        
 
    a_best = chrom(:,:,idx(i)); % top chromosome number i
    disp(a_best)


genetic_offloading_strategy = a_best'.*a; % derive the offloading strategy from the chrromosome

[~,profit_genetic] = feasibility_check(genetic_offloading_strategy);
% [profit_genetic,sol,constraint_check] = get_profit(genetic_offloading_strategy); % get the profit for the offloading strategy using optimization solver 
%     counter =1;
%     end_limit =2;
% while (constraint_check < -1e-2 || profit_genetic == 0) && counter ~= length(a_best) % 
%     for i = 1:end_limit
%     a_best(find(a_best,1,'last'))=0; % remove the last '1' of the chromosome, end_limit times
%     end
%     genetic_offloading_strategy = a_best'.*a;    
%     [profit_genetic,sol,constraint_check] = get_profit(genetic_offloading_strategy);
%     if~mod(counter,3)
%         end_limit = end_limit+1;
%     end
%     
%     counter = counter+1;
% end

profit(i) = profit_genetic;
fprintf('profit_out = %f \n',profit_genetic)
end
profit_genetic = max(profit);  % best of top 5 chromosomes  

end

end
 

function [output]  = reinsert(chrom,selch,objv,objvSel) 
 
[objv_sorted,idx1] = sort(objv,'ascend');
[objvSel_sorted,idx2] = sort(objvSel,'ascend');
for i = 1:length(idx2)
    if(objv_sorted(i)>objvSel_sorted(i))
        output(:,:,i) = chrom(:,:,idx1(i));        
    else
       output(:,:,i) = selch(:,:,idx2(i));
    end
end
    for i = length(idx2)+1:length(idx1)
        output(:,:,i) = chrom(:,:,idx1(i));
    end
    
end
    
    


function [individual] = mutation(individual,mutr)
global M N
individual = reshape(individual,1,[]);
idx = randperm(length(individual),ceil(length(individual)*mutr));
individual(idx)=1-individual(idx);
individual = reshape(individual,M,1);
end

function [cross1 ,cross2]= crossover(parent1,parent2)
global M 
parent1 = reshape(parent1,1,[]);
parent2 = reshape(parent2,1,[]);
total_length = M*1;
split_point= randi([1,M*1]);

cross1(1:split_point) = parent1(1:split_point);
cross1(split_point+1:total_length ) = parent2(split_point+1:total_length );
cross2(1:split_point) = parent2(1:split_point);
cross2(split_point+1:total_length ) = parent1(split_point+1:total_length );
cross1 = reshape(cross1,M,1);
cross2 = reshape(cross2,M,1);
end





function [selch] = rouelttewheelselect(chrom,fitnv,ggap)
[fit,idx] = sort(fitnv,'descend');
S = sum(fit);
[~,~,N] =size(chrom);
new_pop_count = floor(N*ggap);
ret_idx = zeros(1,new_pop_count);
for i=1:new_pop_count
confirm_selection =0;
while(confirm_selection ~= 1)
    
x = rand*S;
j = 1;
P=0;


while(P<=x)
    P = P+fit(j);
    j =j+1;    
end

if(any(ret_idx(:)~=idx(j-1)))
    ret_idx(i) = idx(j-1);
    confirm_selection =1;
end

end
end
selch = chrom(:,:,ret_idx);

end


function [feasible,profit_out] = feasibility_check(a)
% A function which can be used both for checking feasibility and calculating profit for an offloading strategy 
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


if (cvx_status ~= "Infeasible")
    feasible = 1; 
else feasible =0;
end

profit_out = cvx_optval;

end



















