%% Home assignment Derivative pricing Gustaf Sundell, gu0147su-s
clc

n = 12;
sigma = 0.4*ones(n,1);
rho = 0.6;
r = 0.0015;
S0 = 100*ones(n,1);
c = 1/12*ones(n,1);
t = 0;
T = 5;
disc_factor = exp(-r*(T-t));
N = 1000;
K = [80, 100, 120];
SIGMA = get_sigma(sigma,rho,n,T-t);


%% B1

payoffs = zeros(N,1);


prices_1 = zeros(3,1);
for k =1:3
    for i = 1:N

        G = randn(n,1);
        l = log(S0) + (r.*ones(n,1)-diag(SIGMA*SIGMA')/2)*(T-t)+sqrt(T-t)*SIGMA*G;
        R = c'*l;    
        payoffs(i) = max(exp(R)-K(k),0);

    end
    prices_1(k) = disc_factor*mean(payoffs);
end

disp('Lower bound of price (Geometric mean)')
disp('-----------------------------------------------------')
for k = 1:3
   disp('K:')
   disp(K(k))
   disp('Price:')
   disp(prices_1(k))
   disp('-----------------------------------------------------')
end

%% B2

payoffs = zeros(n,N);


prices_2 = zeros(3,1);
for k =1:3
    for i = 1:N

        G = randn(n,1);
        l = log(S0) + (r.*ones(n,1)-diag(SIGMA*SIGMA')/2)*(T-t)+sqrt(T-t)*SIGMA*G;
        S = exp(l);
        payoffs(:,i) = max(S-K(k),0);

    end
    prices_2(k) = c'*disc_factor*mean(payoffs,2);
end
disp('Upper bound of price (average of standard EC:s)')
disp('-----------------------------------------------------')
for k = 1:3
   disp('K:')
   disp(K(k))
   disp('Price:')
   disp(prices_2(k))
   disp('-----------------------------------------------------')
end


%% B3 crude monte carlo


payoffs = zeros(N,1);
prices_3 = zeros(3,1);


for k =1:3
    for i = 1:N

        G = randn(n,1);
        l = log(S0) + (r.*ones(n,1)-diag(SIGMA*SIGMA')/2)*(T-t)+sqrt(T-t)*SIGMA*G;
        S = exp(l);
        payoffs(i) = max(c'*S-K(k),0);

    end    
    prices_3(k) = disc_factor*mean(payoffs);
end

disp('Price of Basket')
disp('-----------------------------------------------------')
for k = 1:3
   disp('K:')
   disp(K(k))
   disp('Price:')
   disp(prices_3(k))
   disp('-----------------------------------------------------')
end



%% B3 - add control variates! 

payoffs = zeros(N,1);
prices_3_control = zeros(3,1);


mean_l = log(S0)+(r*ones(n,1)-diag(SIGMA*SIGMA')/2)*(T-t);
var_l = SIGMA*SIGMA'*(T-t);

mean_R = c'*mean_l;
var_R = c'*var_l*c;

mean_Y = exp(mean_R + 0.5*var_R); % Y = exp(R) log normal... 

Y=zeros(N,1);

for k =1:3
    for i = 1:N

        G = randn(n,1);
        l = log(S0) + (r.*ones(n,1)-diag(SIGMA*SIGMA')/2)*(T-t)+sqrt(T-t)*SIGMA*G;
        S = exp(l);
        Y(i) = exp(c'*l);
        payoffs(i) = max(c'*S-K(k),0);

    end
    b_hat = sum(payoffs.*(Y-mean_Y))/((Y-mean_Y*ones(N,1))'*(Y-mean_Y*ones(N,1)));
    
    prices_3_control(k) = disc_factor*(mean(payoffs) - b_hat.*mean(Y-mean_Y));
end

disp('Price of Basket')
disp('-----------------------------------------------------')
for k = 1:3
   disp('K:')
   disp(K(k))
   disp('Price:')
   disp(prices_3_control(k))
   disp('-----------------------------------------------------')
end

%% Checking conditions
plot(K,prices_1,'-O')
hold on 
plot(K,prices_2,'-O')
hold on 
plot(K,prices_3,'-O')
plot(K,prices_3_control,'-O')
legend({'lower','upper','basket','basket ctrlvariates'})
title('Price at t=0, with varying K on x-axis')

%% B4 
clc
load('HA21_data.mat');
r = -0.001;
T = 5;
t = 0;
n = length(sigma);

sigma = sigma';
S0 = S0';
S1 = S1';

S = diag(sigma);
SIGMA =  chol(S*rho*S,'lower');

N = 10000;

% pr = zeros(N,1);
divs = zeros(N,1);
for i = 1:N
    
   G = randn(n,1);
   l = (r.*ones(n,1)-diag(SIGMA*SIGMA')/2)*(T-t)+sqrt(T-t)*SIGMA*G;
   quotas = exp(l);
   divs(i) = max(c'*quotas -1,0);
%    div = max(c'*quotas -1,0);
%    pr(i) = (1.1*exp(r*(T-t))-1)/div;
    
end
% mean(pr)
pr_hat = (1.1*exp(r*(T-t))-1)/mean(divs)

%% B4 b

t = 1;
phis = zeros(N,1);
NA = 100;

for i = 1:N
    
   G = randn(n,1);
   l = (r.*ones(n,1)-diag(SIGMA*SIGMA')/2)*(T-t)+sqrt(T-t)*SIGMA*G; % a bit unsure if -log(S0) + log(S1) belongs here... 
   quotas = exp(l);
   div = max(c'*quotas -1,0);
   phis(i) = NA*(1+pr_hat*div);
    
end
price_t_1 = exp(-r*(T-t))*mean(phis)

%% funcs

function SIGMA = get_sigma(sigma,rho,n,t)
    if length(sigma)==1
        sigma = sigma.*ones(n,1);
    end
    %assuming rho is scalar for now.
    S = diag(sigma);
    D = rho.*ones(n,n) + diag((1-rho).*ones(n,1));
    SIGMA =  chol(S*D*S,'lower');
end