% dphi function:  derivative of sphi
% Shaobo Han
% 09/24/2015
function output = dphi_f(z, par, opt)
switch opt
    case 'Normal'
        m0 = par(1); v0 = par(2);
        output = -(z./v0).*normpdf(z, m0, sqrt(v0));
end
end