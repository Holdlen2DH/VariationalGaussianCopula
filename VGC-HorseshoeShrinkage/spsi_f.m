%  psi function: small psi, derivative of Psi
% Shaobo Han
% 08/10/2015

function output = spsi_f(x, par, opt)
if nargin<3,
    opt = 'Exp'; % exponential distribution by default
end
if nargin<2,
    par = 0.1; % lambda = 0.1 by default
end
switch opt
    case 'Exp' % exp(lambda)
        output = exppdf(x,1/par); 
    case 'Normal' %
        m0 = par(1); v0 = par(2);
        output =  normpdf(x, m0, sqrt(v0));
    case 'Beta'
        a0 = par(1); b0 = par(2);
        output =  betapdf(x, a0, b0);
end
end