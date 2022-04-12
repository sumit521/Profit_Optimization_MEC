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
% The program execution can take upto 45 min for u=3 (3 iterations) on i7-1170 processor,
% RAM 16GB. More time for higher u. Don't set u=1:1. As the mean() function
% will avarage the single row. Keep u>1, for averaging along columns.
% Have patience!
% Copyright 2021, Sumit Singh, All rights reserved
%%
clear ; close all;

rng(12345)
Fontsize = 12;
global u genetic_fast e omega ohm gamma_C gamma_T Dm Fm B F N M Ln fm_local Tm_max Wm Sm kappa omega pf pt qb qc dataset trainedModel gen best best_feasible best_unpenalized lambda penalty_no 
iter  = 10; % ending UE count 10+iter*delta
start = 2;  % starting UE count 10+start*delta 
delta = 10; % step size for incrementing number of UEs 
load trainedModel.mat 

%% creating random distribution of RRHs 
N = 10; %number of RRHs
cov_rad = 500;
% path loss model 37.6*log(dist) + 148.1
sigma_shadow = 8; %small scale fading model is independently and identically distributed (i.i.d.) Rayleigh fading with zero mean and unit variance
noise_power = -174;
xx0=0; yy0=0; %centre of disk
areaTotal = pi*cov_rad^2; % area of disk

% theta = 2*pi*(rand(N,1));
% rho = cov_rad*sqrt(rand(N,1));
theta = [asin(1/3)*ones(3,1);zeros(4,1);asin(-1/3)*ones(3,1)];
rho = cov_rad;
% [xxRRH, yyRRH] = pol2cart(theta,rho);
% xxRRH = xxRRH + xx0;
% yyRRH = yyRRH + yy0;
xxRRH = [-cov_rad/2,0,cov_rad/2,-3*cov_rad/4,-cov_rad/4,cov_rad/4,3*cov_rad/4,-cov_rad/2,0,cov_rad/2]' + xx0;
yyRRH = [1*cov_rad/2*ones(3,1);zeros(4,1);-1*cov_rad/2*ones(3,1)] + yy0;

for u = 1:10 % number of iterations to average the result
% rng(69); % initialize the seed for resource dependency trends

for q = start:iter


M = 10+delta*q; % number of UEs
%% creating random distribution of UEs
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

%% defining constants
W = zeros(M,3); %task matrix, each row cooresponds to one user's (Fm,Dm and Tm,max)
B = 1e7;    % channel bandwidth 10MHz
Ln = 5e7; % fronthaul capacity for each RRH is 50Mbps
F = 100e9; %MEC server maximum computational capacity 100 GHz
fm_local = 0.7e9;   % local computational capacity 0.7GHz
Dm = unifrnd(50,200,M,1)*1024*8; % output data size
Tm_max = unifrnd(0.6,1.5,M,1)*600e-3; % latency constraint of 600ms
Fm = 1000*Dm; % amountof computation for task Fm
Sm = 0.1; % storage cost 0.1$
kappa = .5;   % impact factor of computation cost
omega = 10;   % impact factor of bandwidth cost
pf = 0.03e-6; % unit price of charge for computation 0.03$/Mega Cycle
pt = 0.3/(1024*1024); %unit price of charge for transmission 0.3$
qb = 0.5e-6;   %unit price of charge for bandwidth 0.5$/Mhz
qc = 0.005e-6; %unit cost for computation 0.005$/mega cycle/s
Wm = [Fm,Dm,Tm_max]; % taks matrix
gamma_T = omega*B*qb*ones(M,1); % communication cost vector
gamma_C = kappa*F*qc*ones(M,1); % computationa cost vector
ohm = (pf*Fm +pt.*Dm+Sm); % revenue vector
dataset = readtable('filtered_dataset.xlsx');
% creating spectrum efficiency matrix  
% e = 5*rand(M,N).*rand(M,N);% spectrum efficiency between each UE and RRH is in range 0 to 5.



%% Calculating the profits 

genetic_fast = 1 ;

for penalty_no = 1:5
   tic
[profit_genetic(u,q,penalty_no),~] = genetic_algo();
time_genetic(u,q,penalty_no) = toc;
end



% for s = 1:2 
%     genetic_fast = s-1; %if genetic_fast == 1, then fast genetic algo is used otherwise normal genetic
%     tic
%     [profit_genetic(u,q,s),~] = genetic_algo();
%     time_genetic(u,q,s)= toc;
% end
% tic
% [profit_sjoora(u,q),S1] = sjoora(); % calculate the profit using sjoora algorithm
% time_sjoora(u,q) = toc;
% 
% 
% tic
% profit_greedy(u,q) = greedy_algo();
% time_greedy(u,q) = toc;
%  tic
% [profit_sbpso(u,q), ~] = SBPSO(); 
% time_sbpso(u,q)= toc;
% tic
% [profit_dsbpso(u,q), ~] = DSBPSO(); 
% time_dsbpso(u,q)= toc;
%%
end
end

%------------------------------------------
% Plotting the results
%------------------------------------------

x = 10+delta*(start:iter); % X axis of the curves, representing number of UEs 

figure;
hold on
legend_text = [];
markers = ['o','+','*','s','d','v','>','h'];
% colors = ['b','c','r','g','k','m'];
colors = [[1 0 0];[0 1 0];[0 0 1];[0.9290 0.6940 0.1250];[0 0 0];[1 1 0]];

for penalty_no = 1:5
    profit_fast_gen(penalty_no,:) = mean(profit_genetic(:,:,penalty_no));
    time_fast_gen(penalty_no,:) = mean(time_genetic(:,:,penalty_no));
    plot(x,profit_fast_gen(penalty_no,start:end),'Color',colors(penalty_no,:),'Marker',markers(penalty_no));
    legend_text{penalty_no} = sprintf('Penalty %d',penalty_no);
end
box on
xlabel("Number of UEs",'FontSize',Fontsize)
ylabel("Profit (in $)",'FontSize',Fontsize)
legend(legend_text,'Location','best');
saveas(gcf,'penaltyComparison','epsc')
saveas(gcf,'penaltyComparison.png')
save PenaltyComparison.mat











%% functions defined for use in script above
%------------------------------------------------------------
% Sticky Binary particle Swarm optimization (SBPSO) Algorithm
%------------------------------------------------------------
function [output_profit,particle] = SBPSO()

global   e N M Wm genetic_fast

disp('Solving using Sticky Binary particle Swarm optimization (SBPSO)')
genetic_fast =1; % for using penalty function


nPop = 15;      % total individuals in population
Varsize = M;    % number of bits in each individual of the population 
T = 30;         % number of iterations
ustkS = 8*T/100;
alpha = 2;
is = 4/M;
ip = alpha*(1-is)/(alpha+1);
ig = (1-is)/(alpha+1);


% create inital offloading strategy matrix 
a = zeros(M,N);
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

%% for finding the number of feasible UEs as starting point
% nones = min(max(dataset.permitted_ues),M);
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




%% create initial population of particles
empty_particle.position = [];
empty_particle.value = 0;
empty_particle.lastposition = [];
empty_particle.prob = rand(1,M);
empty_particle.gbest = [];
empty_particle.pbest = [];
empty_particle.pbestvalue = 0;
empty_particle.gbestvalue = 0;
empty_particle.stk = ones(1,M);




% Create population array
particle = repmat(empty_particle,nPop,1);
% for i = 1:nPop
%     particle(i) = empty_particle;
% end


for j = 1:10
particle(j).position = ismember(1:M,sorted_ues(1:floor((0.5+0.5.*rand)*nones)));

end
for j= 11:nPop
permitted_ues = randsample(1:M,nones,true,rating); %randsample(1:M,nones,true,rating); select permitted UEs based on nones and probability of a UE being selected is proportional to its rating
particle(j).position = ismember(1:M,permitted_ues); % create a chromosome of length M  
end


 [particle(:).lastposition] = (particle(:).position);
 [particle(:).pbest] = particle(:).position;
 [particle(:).gbest] = particle(:).position;


for t=1:T
    %particle(:).gbestvalue = max(particle(:).gbestvalue);
    for i = 1:nPop
        % evaluate particles
        particle(i).value = get_profit_genetic(particle(i).position,a,e);
        if particle(i).value > particle(i).pbestvalue
             particle(i).pbestvalue = particle(i).value;
             particle(i).pbest = particle(i).position;
             if particle(i).value > particle(i).gbestvalue
              [particle(:).gbestvalue] = deal(particle(i).value);  
              [particle(:).gbest] = deal(particle(i).position);
             end
        end
        
        % update gbest,pbest
        
        %update particles
        for j = 1:M
            if particle(i).lastposition(j)~= particle(i).position(j)
                particle(i).stk(j) = 1;
            else
                particle(i).stk(j) = max(particle(i).stk(j)-1/ustkS,0);
            end
            particle(i).prob(j) = is*(1-particle(i).stk(j)) + ip*abs(particle(i).pbest(j)-particle(i).position(j)) + ig*abs(particle(i).gbest(j)-particle(i).position(j));
            
            if rand < particle(i).prob(j)
                particle(i).lastposition(j) = particle(i).position(j);
                particle(i).position(j) = 1-particle(i).position(j);
            end      
            
        end
    end
end

%output = particle;
T = struct2table(particle);
sortedT = sortrows(T,'pbestvalue','descend');
bestof = 3;% for choosing best of n solutions
incr_size = 3;
profit = zeros(1,bestof);
for i =1:bestof
 a_best = sortedT{i,1};
 sbpso_offloading_strategy = a_best'.*a;
 [~,profit(i)] = feasibility_check(sbpso_offloading_strategy);
end

output_profit = max(profit);
while output_profit == -Inf
    for i = bestof:bestof+incr_size
        a_best = sortedT{i,1};
        sbpso_offloading_strategy = a_best'.*a;
        [~,profit(i)] = feasibility_check(sbpso_offloading_strategy);
    end
    output_profit = max(profit);
    bestof = 2*bestof;
    if bestof+incr_size > nPop
        output_profit = 0;
        break
    end
end
disp(profit(:))

end




%------------------------------------------------------------
% Dynamic Sticky Binary particle Swarm optimization (DSBPSO) Algorithm
%------------------------------------------------------------
function [output_profit,particle] = DSBPSO()

global   e N M Wm genetic_fast

disp('Solving using Dynamic Sticky Binary particle Swarm optimization (DSBPSO)')
genetic_fast =1;


nPop = 15;      % total individuals in population
                % number of bits in each individual of the population = M
T = 30;         % number of iterations

ustkSL = 1*T/100;
ustkSU = 8*T/100;
isL = 0/M;
isU = 10/M;


% create inital offloading strategy matrix 
a = zeros(M,N);
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

%% for finding the number of feasible UEs as starting point
% nones = min(max(dataset.permitted_ues),M);
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




%% create initial population of particles
empty_particle.position = [];
empty_particle.value = 0;
empty_particle.lastposition = [];
empty_particle.prob = rand(1,M);
empty_particle.gbest = [];
empty_particle.pbest = [];
empty_particle.pbestvalue = 0;
empty_particle.gbestvalue = 0;
empty_particle.stk = ones(1,M);




% Create population array
particle = repmat(empty_particle,nPop,1);
% for i = 1:nPop
%     particle(i) = empty_particle;
% end


for j = 1:10
particle(j).position = ismember(1:M,sorted_ues(1:floor((0.5+0.5.*rand)*nones)));

end
for j= 11:nPop
permitted_ues = randsample(1:M,nones,true,rating); %randsample(1:M,nones,true,rating); select permitted UEs based on nones and probability of a UE being selected is proportional to its rating
particle(j).position = ismember(1:M,permitted_ues); % create a chromosome of length M  
end


 [particle(:).lastposition] = (particle(:).position);
 [particle(:).pbest] = particle(:).position;
 [particle(:).gbest] = particle(:).position;

% for i = 1:nPop
% particle(i).lastposition = particle(i).position;
% particle(i).pbest = particle(i).position;
% particle(i).gbest = particle(1).position;
% end
for t=1:T
    ustkS = ustkSL + t/T*(ustkSU-ustkSL);
    is = isU - t/T*(isU-isL);
    ig = (1-is)/3;
    ip = 2*ig;
    for i = 1:nPop
        % evaluate particles
        particle(i).value = get_profit_genetic(particle(i).position,a,e);
        if particle(i).value > particle(i).pbestvalue
             particle(i).pbestvalue = particle(i).value;
             particle(i).pbest = particle(i).position;
             if particle(i).value > particle(i).gbestvalue
              [particle(:).gbestvalue] = deal(particle(i).value);  
              [particle(:).gbest] = deal(particle(i).position);
             end
        end
        
        % update gbest,pbest
        
        %update particles
        for j = 1:M
            if particle(i).lastposition(j)~= particle(i).position(j)
                particle(i).stk(j) = 1;
            else
                particle(i).stk(j) = max(particle(i).stk(j)-1/ustkS,0);                
            end
            particle(i).prob(j) = is*(1-particle(i).stk(j)) + ip*abs(particle(i).pbest(j)-particle(i).position(j)) + ig*abs(particle(i).gbest(j)-particle(i).position(j));
            
            if rand < particle(i).prob(j)
                particle(i).lastposition(j) = particle(i).position(j);
                particle(i).position(j) = 1-particle(i).position(j);                
            end            
                
            
        end
    end   
end

%output = particle;
T = struct2table(particle);
sortedT = sortrows(T,'pbestvalue','descend');
bestof = 3;% for choosing best of n solutions
incr_size = 3;
profit = zeros(1,bestof);
for i =1:bestof
 a_best = sortedT{i,1};
 sbpso_offloading_strategy = a_best'.*a;
 [~,profit(i)] = feasibility_check(sbpso_offloading_strategy);
end

output_profit = max(profit);
while output_profit == -Inf
    for i = bestof:2*bestof
        a_best = sortedT{i,1};
        sbpso_offloading_strategy = a_best'.*a;
        [~,profit(i)] = feasibility_check(sbpso_offloading_strategy);
    end
    output_profit = max(profit);
    bestof = 2*bestof;
    if bestof+incr_size > nPop
        output_profit = 0;
        break
    end
end
disp(profit(:))
end









%% SJOORA Algorithm - Refer Algorithm 2 of https://ieeexplore.ieee.org/document/8736324
%---------------------------------------
%   SJOORA Algorithm 
%---------------------------------------

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
%----------------------------
% GREEDY ALGORITHM
%----------------------------
function [profit_out] = greedy_algo()
% function to select the best UEs for offloading
% Refer Algorithm1 in https://ieeexplore.ieee.org/document/6253581 
disp ('Solving using Greedy algorithm')
global   e N M Wm 

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

%% for finding the number of feasible UEs as starting point
% nones = min(max(dataset.permitted_ues),M);
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


%for greedy selection
fractions = [0.95 0.9 0.85 0.8 0.7];

for j = 1:5
best_rated(j,:) = ismember(1:M,sorted_ues(1:ceil(fractions(j)*nones))); % nnz(chrom(:,:,idx(idxmax))))

[~, profit_greedy(j)] =   feasibility_check(best_rated(j,:)'.*a);
end
profit_out = max(profit_greedy);
fprintf('Greedy selection of UEs based on rating. Profit = %f \n',profit_out)

end




%----------------------------------------------------------
% GENETIC ALGORITHM
%----------------------------------------------------------
function [profit_genetic, a_best] = genetic_algo()
% function to select the best UEs for offloading
% Refer Algorithm1 in https://ieeexplore.ieee.org/document/6253581 

global  genetic_fast e N M Wm dataset gen best best_feasible best_unpenalized lambda 

ggap = 1.2; % generation gap
mutr = 0.02; % mutation ratio
a = zeros(M,N);

if genetic_fast == 1
    nind = 30; % value should be greater than 10 and simultaneously also being a multiple of 5. 
    maxgen = 7; % maximum number of generations 
    disp('Solving using Fast Genetic Algorithm')
else
    nind = 15; % if normal genetic, set the the number of individuals in the population low 
    maxgen = 4; % maximum number of generations 
    disp('Solving using Normal Genetic Algorithm')    
end



% defining some vectors for fast genetic algorithm
best = zeros(1,maxgen);
best_feasible = zeros(1,maxgen); 
best_unpenalized = zeros(1,maxgen);
lambda = 25;



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

%% for finding the number of feasible UEs as starting point
% nones = min(max(dataset.permitted_ues),M);
tic
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
toc

for j = 1:10
chrom(:,:,j) = ismember(1:M,sorted_ues(1:floor((0.5+0.5.*rand)*nones))); % 

end

%% create initial population

for j= 11:nind
permitted_ues = randsample(1:M,nones,true,rating); %randsample(1:M,nones,true,rating); select permitted UEs based on nones and probability of a UE being selected is proportional to its rating
chrom(:,:,j) = ismember(1:M,permitted_ues); % create a chromosome of length M  
end


gen_obs_mut = 2;
gen_obs_break = 3;
gen =1;


while(gen < maxgen) % while the maximum number of generations not reached
    
    
    %---------------------------------------------------
% parameters update for adaptive penalty functions  
Nf =3; % previous generations to consider
beta1 = 3;
beta2 = 2;
if mod(gen,Nf) == 0
    if sum(best_feasible(gen-Nf+1:gen)) == 0
        lambda = lambda*beta1;
    elseif sum(best_feasible(gen-Nf+1:gen)) == Nf
        lambda = lambda/beta2;
    end
end
%----------------------------------------------------  
  
    %for genetic algorithms parameters control on fly
    if gen > gen_obs_break
        
        if mod(gen,gen_obs_mut) == 0 && std(info_gen(gen-gen_obs_mut:gen-1,2)) < 5
            mutr = max(mutr + 0.01,0.01);% to avoid setting it to zero
            fprintf('Changing mutation rate. current rate = %f \n',mutr)
        end
        
        if mod(gen,gen_obs_break) == 0 && std(info_gen(gen-gen_obs_break:gen-1,2)) < 5
            fprintf('Breaking out of the evolution. std deviation of best fitness over last %d generations = %f \n',gen_obs_break,std(info_gen(gen-gen_obs_break:gen-1,2)))
            break
        end
    end
    

    for i = 1:nind
      objv(i) = get_profit_genetic(chrom(:,:,i),a,e); % calculate the profit using the function. Create objctive vector of profits 
    end       
        info_gen(gen,:) = [mean(objv), max(objv)];
        fitnv = objv./norm(objv); % create normalized fitness vector
       new_nind = floor(ggap*nind); % new individuals to be created for next generation
       
%       parentch =rouelttewheelselect(chrom,fitnv,ggap); % roulette wheel parent selection
       
       parentch = tournamentselection(chrom,fitnv,new_nind); % tournament based parent selection       
        
        i=1;
        while(i<new_nind )
            [chldch(:,:,i),chldch(:,:,i+1)] = crossover(parentch(:,:,i),parentch(:,:,i+1)); % crossover of selected chromosomes
            i = i+2;
        end
        
        for i=1:new_nind
            mutatedch(:,:,i) = mutation(chldch(:,:,i),mutr); % mutation of selected chromosomes s per mutation ratio
        end
        
        for i=1:new_nind
            objvmutatedch(i) = get_profit_genetic(mutatedch(:,:,i),a,e); % again created an objective vector of profits for selected chromosomes          
        end
%        chrom = reinsert(chrom,mutatedch,objv,objvmutatedch); % reinsert the new individuals in the population by slecting best among old and new such that toatl count remains nind   
        chrom = lambda_mu_selection(mutatedch,objvmutatedch,nind);
        gen = gen+1;  
    
end
% figure;
% subplot(1,2,1)
% plot(info_gen(:,1))
% xlabel("Generation",'FontSize',12)
% ylabel(" Profit (in $)",'FontSize',12)
% title('mean over generations')
% subplot(1,2,2)
% plot(info_gen(:,2))
% xlabel("Generation",'FontSize',12)
% ylabel("Profit (in $)",'FontSize',12)
% title('Approximate maximum fitness over generations')

    for i = 1:nind
        objv(i) = get_profit_genetic(chrom(:,:,i),a,e);% calculate the profits for the final population after nind generations
    end
  
    [objv_sorted,idx] = sort(objv,'descend'); % sort the profits in descending order
 
if genetic_fast ~= 1    
    profit_genetic = objv_sorted(1); % if it is normal genetic algo, set the output profit as top of sorted objv 
    a_best = idx(1); % best chromosome corresponding to chromosome with highest profit
    else %condition for fast genetic algorithm
%%  for fast genetic algorithm

bestof = 3;% for choosing best of n solutions
incr_size = 3;
profit = zeros(1,bestof);
for i =1:bestof
 a_best = chrom(:,:,idx(i)); % top chromosome number i;
 genetic_offloading_strategy = a_best'.*a;
 [~,profit(i)] = feasibility_check(genetic_offloading_strategy);
end

profit_genetic = max(profit);
while profit_genetic == -Inf
    for i = bestof:2*bestof
        a_best = chrom(:,:,idx(i)); % top chromosome number i;
        genetic_offloading_strategy = a_best'.*a;
        [~,profit(i)] = feasibility_check(genetic_offloading_strategy);
    end    
    profit_genetic = max(profit);
    bestof = 2*bestof;
    if bestof+incr_size > nind
        profit_genetic = 0;
        break
    end
end
disp(profit(:))


end
fprintf('Final profit from genetic algorithm = %f \n',profit_genetic)
toc

end


function [output] = lambda_mu_selection(population,objpopulation,nind)
[~,idx] = sort(objpopulation,'descend');
output = population(:,:,idx(1:nind));
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
individual = reshape(individual,1,[]);

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

function [selch]= tournamentselection(chrom,fitnv,new_nind)
[~,idx]=sort(fitnv,'descend');
i=1;
nind = length(idx);
k = 0.2*nind; % 20% of population to be selected for toutnament. Reduce the value for higher population size, i.e. higher nind.

while  i <= floor(new_nind)
    tournament(:,:,i) = randperm(nind,k);  
    [~,loc] = ismember(tournament(:,:,i),idx);
    winner = min(loc);
    selch(:,:,i) = chrom(:,:,winner);
    i = i+1;
end

end


function [feasible,profit_out] = feasibility_check(a)
% A function which can be used both for checking feasibility and calculating profit for an offloading strategy 
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



















