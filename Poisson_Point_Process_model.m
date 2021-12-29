% https://www.researchgate.net/post/How_can_I_model_a_simple_communication_channel_using_MATLAB
% https://www.researchgate.net/post/Who-knows-how-to-model-a-simple-communication-channel-using-MATLAB
% https://hpaulkeeler.com/simulating-a-poisson-point-process-on-a-disk/
close;
N = 10;

M = 50;
cov_rad = 500;
% path loss model 37.6*log(dist) + 148.1
sigma_shadow = 8; %small scale fading model is independently and identically distributed (i.i.d.) Rayleigh fading with zero mean and unit variance
noise_power = -174;
xx0=0; yy0=0; %centre of disk
areaTotal = pi*cov_rad^2; %area of disk

lambda = 10;
% numbpoints = poissrnd(areaTotal*lambda);
numbPoints = N;
% theta = 2*pi*(rand(numbPoints,1));
% rho = cov_rad*sqrt(rand(numbPoints,1));
theta = [asin(1/3)*ones(3,1);zeros(4,1);asin(-1/3)*ones(3,1)];
rho = cov_rad;
% [xxRRH, yyRRH] = pol2cart(theta,rho);
% xxRRH = xxRRH + xx0;
% yyRRH = yyRRH + yy0;
xxRRH = [-cov_rad/2,0,cov_rad/2,-3*cov_rad/4,-cov_rad/4,cov_rad/4,3*cov_rad/4,-cov_rad/2,0,cov_rad/2]' + xx0;
yyRRH = [1*cov_rad/2*ones(3,1);zeros(4,1);-1*cov_rad/2*ones(3,1)] + yy0;
scatter(xxRRH,yyRRH, 'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5);
xlabel('x');ylabel('y');

axis square;
hold on
thetaUE = 2*pi*(rand(M,1));
rhoUE = cov_rad*sqrt(rand(M,1));

[xxUE, yyUE] = pol2cart(thetaUE,rhoUE);
xxUE = xxUE + xx0;
yyUE = yyUE + yy0;
scatter(xxUE,yyUE,'d')
legend('RRH','UE')
%title('Location of RRH and UEs ')
hold off
exportfig(gcf,'RRHlocation.eps','Color','rgb','LineWidth',2)
saveas(gcf,'RRHlocation.png')
saveas(gcf,'try','epsc')
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
        
var(max(e,[],2))
channel = sqrt(variance/2)*(randn + randn);





