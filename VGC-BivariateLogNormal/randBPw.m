% Generate a P*D BP weight matrix
% (every column is a probability vector which sums up to one)
% P: number of variables
% D: number of BP basis functions
% a: shape parameter; b: scale parameter

% Shaobo Han
% 08/08/2015
function w = randBPw(P, D, a, b)
if nargin<4
    a = 1; b = 1; % by default
end
if D==1
    w = ones(P,1);
else
    w0 = gamrnd(a,b,P,D);
    w = diag(1./sum(w0,2))*w0;
end
end