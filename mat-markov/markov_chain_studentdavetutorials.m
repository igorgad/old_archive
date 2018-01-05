

% written by StudentDave
%for licensing and usage questions
%email scienceguy5000 at gmail. com

%% Welcome to the BAYESIAN DOJO! Here, we'll learn about markov chains

% our main examples will be of Ergodic regular markov chains
% these type of chains converge to a steady-state, and have some nice
% properties for rapid calculation of this steady state.


%% non-regular
%first, an example of a non-regular, but ergodic markov chain, which
%doesn't converge

P = [ 0 1 ;1 0]

t_all = [];
i_all = [];
figure(1)
clf
for i = 1:100
    t = P^i
    t_all = [t_all t(:)];
    i_all = [i_all  ones(4,1)*i];
    subplot(211)
    draw_states(t,i)
    subplot(212)
    plot(i_all',t_all','.-')
    xlabel('discrete time steps')
    ylabel('probability')
    title('evolution of transition probs. for each element')
    pause
end



%% Has zero elements in initial transition matrix (TM)), but not always, so
% thus still regular
P = [ 1/2 1/2 ;1 0]


t_all = [];
i_all = [];
figure(1)
clf
for i = 1:100
    t = P^i
    t_all = [t_all t(:)];
    i_all = [i_all  ones(4,1)*i];
    subplot(211)
    draw_states(t,i)
    subplot(212)
    plot(i_all',t_all','.-')
    xlabel('discrete time steps')
    ylabel('probability')
    title('evolution of transition probs. for each element')
    pause
end

%% Training: level one! Learn to read your opponent
% Our Bayesian Ninja has just learned about the sudden reappearance of an
% old enemy of the Bayesian Clan, The frequentisian Ninja Clan!
% The bayesian ninja must now train to fight!
% One critical skill is to know your opponent, know their tendencies, and
% learn their patterns!  The old master ninja is the last of the Bayesians who have fought the
% Frequentisian Ninja, and he describes the different close-combat fighting styles
% taught within their devious clan in terms of a markov process.

% In the first lesson, he describes a three state fighting style, comprised of
% 1)punch (red) and 2)kick (yellow), 3) and flying falcon punch (blue).
% here we look at how, overall, the Frequentisian Ninja's will fight, given
% the probabilities of how they mix up their punching, kicking, and
% "special attack"
%

% E honda style (Likes to punch)-------------------------------

a= .9
b = .3
c = .2
P = [ a (1-a)/2 (1-a)/2;  (1-b)/2 b (1-b)/2; (1-c)/2 (1-c)/2 c ;]


% chun Li Style (Likes to kick)-------------------------------

a= .1
b = .9
c = .3
P = [ a (1-a)/2 (1-a)/2;  (1-b)/2 b (1-b)/2; (1-c)/2 (1-c)/2 c ;]

% Captain Falcon Style (Obvious? :) -------------------------------

a= .1
b = .2
c = .9
P = [ a (1-a)/2 (1-a)/2;  (1-b)/2 b (1-b)/2; (1-c)/2 (1-c)/2 c ;]

% MASTER Frequentisian Ninja (Perfectly equal skill in each

a = .33333
b = .333333
c = .333333
P = [ a (1-a)/2 (1-a)/2;  (1-b)/2 b (1-b)/2; (1-c)/2 (1-c)/2 c ;]

t_all = [];
i_all = [];
figure(1)
clf
for i = 1:100
    t = P^i
    t_all = [t_all t(:)];
    i_all = [i_all  ones(size(t_all,1),1)*i];
    subplot(211)
    draw_states3(t(:),i)
    subplot(212)
    plot(i_all',t_all','.-')
    xlabel('discrete time steps')
    ylabel('probability')
    title('evolution of transition probs. for each element')
    axis([0 max(max(i_all)) min(min(t_all))-.5 max(max(t_all))+.5])
    pause
end



%% Training: level two! What's it take to beat the MASTER Frequentisian Ninja?
% here, we show how, even with just a 1st order markov chain, you can still
% simulation systems that depend on past events.
% If the Master ninja can land a 3 hit combo in the specific order of
% punch, kick, falcon punch...the bayesian ninja will get KO'd!
% Simulate this as a markov process and see how the Bayesian will perform
% given his ability, b, to interrupt the punches, and starting state u,
% after T number of attacks (time steps)


a= .5
b = .7 %interrupt probability
P = [1-a a 0 0; b 0 1-b 0; b 0 0 1-b; 1 0 0 0]

t_all = [];
i_all = [];
figure(1)
clf

for i = 1:100
    t = P^i
    t_all = [t_all t(:)];
    i_all = [i_all  ones(size(t_all,1),1)*i];
    subplot(211)
    draw_states4(t,i)
    subplot(212)
    plot(i_all',t_all','.-')
    xlabel(['Time steps = ', num2str(i)])
    ylabel('probability')
    title('evolution of transition probs. for each element')
    pause
end


%how to solve with eigen math :)--------------------

% get eigen decomposition 
[Evector,Evalue] = eig(P')
%get values out of matrix
values     = diag(Evalue);
%find the unitary Evalue, won't be exact, so use generic tool for finding
%closest element to N. There will only be one for regular transition matrix

N = 1
[min_v,coln] = min(abs(values-N))
%grab the corresponding vector, and normalize, this is your stationary
%distribution!
Evector = Evector(:,coln)

fixed_row_vector = (Evector/sum(Evector))'


