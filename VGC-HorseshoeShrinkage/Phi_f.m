% Phi function: mapping from (-infty, infty) to unit closed interval [0,1]
% Shaobo Han
% 09/17/2015
function u = Phi_f(z, par, opt)
switch opt
    case 'Normal'
        m0 = par(1); v0 = par(2);
        u = normcdf(z, m0, sqrt(v0));
end
end
