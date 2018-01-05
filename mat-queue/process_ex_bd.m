% process_ex_bd.m sets up parameters for a simple birth/death process

Tmax   = 5;                             % maximum time
beta   = 1.2;                           % individuals' birth rate
delta  = 1;                             % individuals' death rate
p      = beta/(beta+delta);             % probability of a birth at jump time
n      = 20;                            % maximum population size

S      = 0:n;                           % state space
mu     = zeros(1,n+1);                  % initial distribution
mu(min(5,n)) = 1;
lambda =(beta+delta)*(0:n);             % sojourn parameters
Q=diag(p*ones(1,n),1)+diag((1-p)*ones(1,n),-1); % transitions up or down
Q(n+1,n+1)=p;                           % limit population size at n
Q(1,1)=1;                               % make sure 0 is absorbing
Q(1,2)=0;
lambda(1)=1;                            % avoid problems with simulation