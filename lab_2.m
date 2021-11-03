%% lab 2 Gustaf Sundell & Fredrik Lindstedt


%% Basket call option 2

K=100;
T=4;
r = 0.0015;
St = 100;
sigma = 0.4;
rho=0.6;
n=12;
N=10000;
t=0;
clc

S = diag(sigma);
D = rho.*ones(n,n) + diag((1-rho).*ones(n,1));
SIGMA = S*D*S*T;
SIGMA = chol(SIGMA,'lower'); 

ST_basket = zeros(N,n);
payoffs_basket = zeros(N,1);
% simulating payoffs N times
for i =1:N
    G = randn(n,1);
    ST_basket(i,:) = St.*exp((r.*ones(n,1)-diag(SIGMA*SIGMA')/2)*(T-t)+SIGMA*G*sqrt(T-t));
    payoffs_basket(i) = max(mean(ST_basket(i,:)-K), 0);
end

disc_factor = exp(-r*(T-t));
disp('Basket option price:')
Pi_hat = mean(payoffs_basket).*disc_factor
disp('Standard deviation:')
dev = std(payoffs_basket*disc_factor)/sqrt(N)

%% Heston model 2
clc

K = 80;
St = 90;
T = 1;
r = 0.01;
kappa = 10;
theta = 0.16;
sigma = 0.1;
rho = -0.8;
V0 = 0.16;
C = 1;
p = 1;
par = [V0,kappa,theta,sigma,rho];

disp('opt_price:')
P=opt_price('Heston',par,C,p,K,St,r,T)

steps = 500;
N = 10000;
h = T/steps;

d = 4*theta*kappa/(sigma^2);
lambda = 4*kappa*exp(-h*kappa)/(sigma^2*(1-exp(-h*kappa)));
C = sigma^2*(1-exp(-h*kappa))/(4*kappa);

sim_V = zeros(N+1,steps);
sim_S = zeros(N+1,steps);

sim_V(:,1) = V0;
sim_S(:,1) = St;
G = randn(N+1,steps);

for i = 2:steps
    sim_V(:,i) = C*ncx2rnd(d,sim_V(:,i-1).*lambda);
    sim_S(:,i) = sim_S(:,i-1).*exp(h.*(r-rho.*kappa.*theta./sigma + 0.5.*(sim_V(:,i) + sim_V(:,i-1)).*(kappa.*rho./sigma - 0.5)) + (rho./sigma).*(sim_V(:,i) - sim_V(:,i-1)) + sqrt(h.*0.5.*(sim_V(:,i) + sim_V(:,i-1)).*((1-rho).^2)).*G(:,i));
end

payoffs_heston = max(sim_S(:,end) - K, 0);

disc_factor = exp(-r*(T-t));
disp('European call option price using Heston method:')
Pi_hat = mean(payoffs_heston).*disc_factor
disp('Standard deviation:')
dev = std(payoffs_heston*disc_factor)./sqrt(N)

%sim_S(1:10,end)







%% hej

K=80;
T=1;
r = 0.01;
St = 90;
sigma = 0.4;
rho=0.6;
n=12;
N=10000;
t=0;
clc

sim_all(r,T,t,N,K,St,sigma,rho,n)
if n == 1
    disp('price by bls')
    bls=blsprice(St, K, r, T, sigma)
end

sim_all(r,T,t,N,K,St,sigma,rho,n,'antithetic')
if n == 1
    disp('price by bls')
    bls=blsprice(St, K, r, T, sigma)
end
sim_all(r,T,t,N,K,St,sigma,rho,n,'control') 
if n == 1
    disp('price by bls')
    bls=blsprice(St, K, r, T, sigma)
end

%% basket option


K=80;
T=1;
r = 0.0015;
St = 100;
sigma = 0.4;
rho=0.6;
n=12;
N=10000;
t=0;
clc

sim_all(r,T,t,N,K,St,sigma,rho,n)
if n == 1
    disp('price by bls')
    bls=blsprice(St, K, r, T, sigma)
end

sim_all(r,T,t,N,K,St,sigma,rho,n,'antithetic')
if n == 1
    disp('price by bls')
    bls=blsprice(St, K, r, T, sigma)
end
sim_all(r,T,t,N,K,St,sigma,rho,n,'control') 
if n == 1
    disp('price by bls')
    bls=blsprice(St, K, r, T, sigma)
end


%% Heston model
clc
kappa = 10;
theta = 0.16;
sigma = 0.1;
rho = -0.8;
V0 = 0.16;
C = 1;
p = 1;
Tau = T;
par = [V0,kappa,theta,sigma,rho];


disp('opt_price:')
P=opt_price('Heston',par,C,p,K,St,r,Tau)

steps = 1000;
N = 1000;
h = T/steps;

d = 4*theta*kappa/(sigma^2);
lambda = 4*kappa*exp(-h*kappa)/(sigma^2*(1-exp(-h*kappa)));
C = sigma^2*(1-exp(-h*kappa))/(4*kappa);

sim_V = zeros(N+1,steps);
sim_S = zeros(N+1,steps);

sim_V(:,1) = V0;
sim_S(:,1) = log(St);
G = randn(N+1,steps);

for i =2:steps
    sim_V(:,i) = C*ncx2rnd(d,sim_V(:,i-1).*lambda);
    
    tmp = sim_S(:,i-1) + h*(r-rho*kappa*theta/sigma); 
    tmp =  tmp + h*0.5*(sim_V(:,i)+sim_V(:,i-1))*(kappa*rho/sigma - 0.5);
    tmp = tmp + (rho/sigma)*(sim_V(:,i)-sim_V(:,i-1));
    tmp = tmp + sqrt(h*0.5*(sim_V(:,i)+sim_V(:,i-1))*(1-rho^2)).*G(:,i);
    
    sim_S(:,i) = tmp;
    
end

ST = exp(sim_S(:,end));

price_sim_heston = exp(-r*(T-t))*mean(EC_payoff(ST,K,1))

std_heston = std(EC_payoff(ST,K,1))/sqrt(N)

%% forts

S = exp(sim_S);
m = max(S')';
barr = exp(-r*(T-t))*mean(EC_payoff(ST,K,1).*(m<100))


%% last task
clc
St = 50;
B = 100;
K = 40;

tot_steps=1000;
N = 2000;
S = zeros(N,tot_steps);
delta = T/tot_steps;
S(:,1) = St;
payoffs = zeros(N,1);

for sim = 1:N
    for step = 2:tot_steps
     s = S(sim,step-1);   
     G = randn();   
     S(sim,step) = s + r*s*delta + sigma*s*sqrt(delta)*G + sigma^2*s*0.5*h*(G^2-1);
     
    end
    payoffs(sim) = barrier_payoff(S(sim,:),K,B);
end
disc_factor = exp(-r*(T-t));    
disp('price for barrier option')
price_barrier = disc_factor*mean(payoffs)
dev_barrier = std(payoffs)/sqrt(N)

disp('price by bls')
bls=blsprice(St, K, r, T, sigma)

%% functions

function phi = barrier_payoff(S_trad,K,B)
    if max(S_trad)>= B
        phi = 0;
    else
        phi = max(S_trad(end)-K,0);
    end

end

function SIGMA = get_cov_mat(sigma,rho,n,t)
    if length(sigma)==1
        sigma = sigma.*ones(n,1);
    end
    %assuming rho is scalar for now.
    S = diag(sigma);
    D = rho.*ones(n,n) + diag((1-rho).*ones(n,1));
    SIGMA =  S*D*S*t;


end

function PHI = EC_payoff(S,K,n)
    if n == 1
        PHI = max(S-K,0);
    else
        PHI = max(mean(S)-K,0);
    end
end

function PHI = simulate_ST(St,r,n,Sigma,T,t,method,K)
    if strcmp(method,'standard')
        G = randn(n,1);
%         St.*exp((r.*ones(n,1)-diag(Sigma*Sigma')/2)*(T-t)+Sigma*G*sqrt(T-t))
        PHI = EC_payoff(St.*exp((r.*ones(n,1)-diag(Sigma*Sigma')/2)*(T-t)+Sigma*G*sqrt(T-t)),K,n);
        
    elseif strcmp(method,'antithetic')
        u = rand(n,1);
        G1 = norminv(u);
        PHI1 = EC_payoff(St.*exp((r.*ones(n,1)-diag(Sigma*Sigma')/2)*(T-t)+Sigma*G1*sqrt(T-t)),K,n);
%         PHI1
        G2 = norminv(1-u);
        PHI2 = EC_payoff(St.*exp((r.*ones(n,1)-diag(Sigma*Sigma')/2)*(T-t)+Sigma*G2*sqrt(T-t)),K,n);
        
        PHI = (PHI1+PHI2)/2;
%         PHI

    elseif strcmp(method,'control')
        G = randn(n,1);
        PHI = St.*exp((r.*ones(n,1)-diag(Sigma*Sigma')/2)*(T-t)+Sigma*G*sqrt(T-t));
        
    else
    end

end

function [Pi_hat,dev] = est_Pi(payoffs,r,T,t)
    disc_factor = exp(-r*(T-t));
    Pi_hat = mean(payoffs,1).*disc_factor;
    dev = std(payoffs);
end

function [Pi_hat,dev] =control_eval(ST,St,r,T,t,K,N,n)
    disc_factor = exp(-r*(T-t));
    Y = ST.*disc_factor;
    
    payoffs = EC_payoff(ST,K,n);
    
  
    b_hat = sum(payoffs.*(Y-St),1)/sum((Y-St).^2,1);
    Pi_hat = disc_factor*mean(payoffs)-(b_hat/N).*sum(mean(Y,2)-St,1);
    cv = cov(payoffs,Y);
    var_control = var(payoffs)-cv(1,2)^2/(N*var(Y));
    dev = sqrt(var_control);
end

function sim_all(r,T,t,N,K,St,sigma,rho,n,method)
    if nargin <10
        method = 'standard';
    end
    
    if n == 1
        SIGMA = sigma;
    else
        SIGMA = get_cov_mat(sigma,rho,n,T-t)
    end
    
    if strcmp(method,'control')
       ST = zeros(N,n);
       for i =1:N
          ST(i,:) = simulate_ST(St,r,n,SIGMA,T,t,method,K)'; 
       end
       disp(strcat('price by simulating using ',method))
       [Pi_hat,dev]=control_eval(ST,St,r,T,t,K,N,n);
       Pi_hat
       disp('and std.dev:')
       dev
       disp('if div by sqrt(N):')
       dev/sqrt(N)
       
    else

        PHI = zeros(N,1);
        for i =1:N
            PHI(i) = simulate_ST(St,r,n,SIGMA,T,t,method,K)';
        end
    %     PHI

        disp(strcat('price by simulation using ',' ',method))
        [Pi_hat,dev] = est_Pi(PHI,r,T,t);
        Pi_hat
       disp('and std.dev:')
       dev
       disp('if div by sqrt(N):')
       dev/sqrt(N)
    end

end