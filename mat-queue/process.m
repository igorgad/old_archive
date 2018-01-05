% process.m simulates a Markov process on a state space S, using initial
% distribution mu, sojourn parameters lambda, and jump matrix Q.
clear;
close all;

Tmax   = 100;                             % maximum time
beta   = 1.1;                           % individuals' birth rate
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

clear T;
clear x;             % clear out previous values

T(1) = 0;            % start times at 0
x(1) = rando(mu);    % generate first x value (time 0, not time 1)
i    = 1;

while T(i) < Tmax,
  T(i+1) = T(i)  - log(rand)/lambda(x(i));  % generate exponential rv
                                % occasionally causes errors; ignore them
  x(i+1) = rando(Q(x(i),:));    % use Q to make state transitions
  i=i+1;
end

cla;
hold on
% 
% for i=1:(length(T)-1),
%   plot([T(i) T(i+1)], [S(x(i)) S(x(i))]);
% end

stairs(T,S(x));                 % sometime this is more appropriate
axis([0 Tmax 0 (length(mu)+1)]);
xlabel('Time');
