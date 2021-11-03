%% labb 1 Derivatie pricing, Gustaf and Fredrik
sigma = 0.2;
T= 1;
n = 3;
delta = T/n;
r = 0.05;

K = 100;
S0 = 90;

u = exp(sigma*sqrt(delta));
d = 1/u;

qu = (exp(r*delta)-d)/(u-d);
qd=1-qu;

binmat = zeros(n+1,n+1);


for j=1:n+1
    for i = 1:j
       binmat(i,j) = S0*u^(j-i)*d^(i-1); 
    end
end
%% RNVF
price = 0;
for k=0:n
    price = price + nchoosek(n,k)*qu^(n-k)*qd^(k)*payoff_EC(K,binmat(k+1,end));
end

price = price*exp(-r*T);
price

%% using nodes
pricemat= zeros(n+1,n+1);
pricemat(:,end) = payoff_EC(K,binmat(:,end));

for k=n:-1:1
    for j = 1:k
       pricemat(j,k) = exp(-r*delta)*(qu*pricemat(j,k+1) + qd*pricemat(j+1,k+1)); 
        
    end
end
C = pricemat;
C
%% american put
pricemat = zeros(n+1,n+1);
pricemat(:,end) = payoff_AP(K,binmat(:,end));

for k=n:-1:1
    for j = 1:k
       pricemat(j,k) = max(exp(-r*delta)*(qu*pricemat(j,k+1) + qd*pricemat(j+1,k+1)),payoff_AP(K,binmat(j,k))); 
        
    end
end
AP = pricemat;
AP
%% Forward
pricemat= zeros(n+1,n+1);
pricemat(:,end) = binmat(:,end)-K;


for k=n:-1:1
    for j = 1:k
       pricemat(j,k) = exp(-r*delta)*(qu*pricemat(j,k+1) + qd*pricemat(j+1,k+1)); 
        
    end
end
F = pricemat;
F

%% condition check

bool = binmat-K <= C-AP & C-AP <= F

%% 3.3 convergence

N = 100;
allprices = zeros(1,N);
for n = 1:N
    delta = T/n;
    u = exp(sigma*sqrt(delta));
    d = 1/u;
    qu = (exp(r*delta)-d)/(u-d);
    qd=1-qu;
    
    % defining binmat
    binmat = zeros(n,n+1);
    for i=1:n+1
        for j = i:n+1
            binmat(i,j) = S0 * u^(j-i) * d^(i-1);
        end
    end
    % RNVF
    price = 0;
    for k=0:n

        price = price + nchoosek(n,k)*qu^(n-k)*qd^k*payoff_EC(K,binmat(k+1,end));
    end
    price = price*exp(-r*T);
    allprices(i) = price;
end
allprices;

d1 = (log(S0/K)+(r+sigma^2/2)*T)/(sigma*sqrt(T));
d2 = d1-sigma*sqrt(T);

BS_price = S0*normpdf(d1)-K*exp(-r*T)*normpdf(d2);

plot(allprices)
hold on
yline(BS_price)

%% functions

function phi = payoff_EC(K,ST)
    phi = max(ST-K,0);
end

function phi = payoff_AP(K,St)
    phi = max(K-St,0);
end