% phi function: small phi, derivative of Phi
% Shaobo Han
% 08/07/2015
function output = sphi_f(z, par, opt)
switch opt
    case 'Normal'
        m0 = par(1); v0 = par(2);
        output = normpdf(z, m0, sqrt(v0));
end
end
