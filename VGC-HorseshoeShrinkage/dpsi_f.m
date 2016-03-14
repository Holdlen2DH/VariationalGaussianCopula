%  dpsi function: derivative of small psi
% Shaobo Han
% 08/10/2015

function output = dpsi_f(x, par, opt)
if nargin<3,
    opt = 'Exp'; % exponential distribution by default
end
if nargin<2,
    par = 0.1; % lambda = 0.1 by default
end
switch opt
    case 'Exp' % exp(lambda)
        output = -par.*exppdf(x,1/par);
    case 'Normal' %
        m0 = par(1); % mean
        v0 = par(2); % variance
        output =  (-x./v0).*normpdf(x,m0,sqrt(v0));
    case 'Beta'
        a0 = par(1); b0 = par(2);
        output =(a0+b0-1)*(betapdf(x, a0-1, b0)- betapdf(x, a0, b0-1));
end
end