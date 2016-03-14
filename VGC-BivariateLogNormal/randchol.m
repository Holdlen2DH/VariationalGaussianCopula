% Generate a P*P lower triangular Cholesky factor matrix
% Shaobo Han
% 08/07/2015
function C = randchol(P)
a = rand(P, P);
aTa = a'*a;
[V, D] = eig(aTa);
aTa_PD = V*(D+1e-5)*V';
C=chol(aTa_PD, 'lower');
end