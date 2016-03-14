% Invserse Psi function:
% Shaobo Han
% 08/10/2015

function output = iPsi_f(p, par, opt)
if nargin<3,
    opt = 'Exp'; % exponential distribution by default
end
if nargin<2,
    par = 0.1; % lambda = 0.1 by default
end

p(p>1-1e-16) = 1-1e-16; % numerical stability 

switch opt
    case 'Exp' % exp(lambda)
        %    output = -log(1-p)./par;
        output = expinv(p,1/par); 
    case 'Normal' %
        m0 = par(1); v0 = par(2);
        output =  norminv(p,m0,sqrt(v0));
    case 'Beta'
        a0 = par(1); b0 = par(2);
        output =  betainv(p, a0, b0);
end
end