% graph_matrix.m produces a graphical display of a transition matrix
% a, b are the lowest and highest states of the chain
% P is the transition matrix

function [void] = graph_matrix(a, b, P)

x=(a-0.5):1:(b+0.5);                               % column numbers

M = length(P(1,:));                                % number of states

if M == (b-a+1),                                   % P is the right size
  pcolor(x, x, [[P zeros(M,1)]' zeros(M+1,1)]');   % display matrix P
  axis ij;          % number the axes as for a matrix, not a regular graph
end