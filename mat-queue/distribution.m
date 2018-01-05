% distribution.m computes and graphs the distribution of a Markov process
% over time

Tmax   = 5;                             % maximum time
beta   = 23.8;                           % individuals' birth rate
delta  = 24;                             % individuals' death rate
p      = beta/(beta+delta);             % probability of a birth at jump time
n      = 5;                            % maximum population size

S      = 0:n;                           % state space
mu     = zeros(1,n+1);                  % initial distribution
mu(min(5,n)) = 1;
lambda =(beta+delta)*(0:n);             % sojourn parameters
Q=diag(p*ones(1,n),1)+diag((1-p)*ones(1,n),-1); % transitions up or down
Q(n+1,n+1)=p;                           % limit population size at n
Q(1,1)=1;                               % make sure 0 is absorbing
Q(1,2)=0;
lambda(1)=1;                            % avoid problems with simulation

% Simply modify the initial distribution, lambda, jump matrix Q, and
% state space S

N     = length(S);                        % number of states
step  = Tmax/100;                         % standard time increment
T     = 0:step:Tmax;                      % vector of times

A     = diag(lambda)*(Q-eye(size(Q)));    % generator
Pstep = expm(A*step);                     % this is P(step); matrix exponential
m     = mu;                               % distribution at each step
dist  = mu;                               % each row of dist is distn at a time

for t=1:(length(T)-1),
  m    = m*Pstep;                         % step forward one time step
  dist = [dist; m];                       % add another row to m
end

for v=1:N,
  subplot(N,1,N-v+1);                     % one plot for each state
  plot(T,dist(:,v));                      % plot probability over time
  axis([0 Tmax 0 1]);
  ylabel(['P(X_{t} = ' num2str(S(v)) ')']);
end

subplot(N,1,1);
title('Probabilities of being in states 1, 2, 3, ... over time');

% Now compute the invariant/limiting distribution

eta = invariant(expm(A));

eta2 = invariant(Q)*diag(1./lambda);
eta2 = eta2/sum(eta2);                  % normalize


