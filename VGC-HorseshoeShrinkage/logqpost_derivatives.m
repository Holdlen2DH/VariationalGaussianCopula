% Derivatives of Log likelihood and Log prior
% Model dependent term
% Shaobo Han
% 08/10/2015

function g = logqpost_derivatives(x, model, par)
switch model
    case 'Horseshoe'
        y = par.y;
        g = zeros(2,1);
        g(1) = -2/x(1)+y^2/2/(x(1)^2)+x(2)/(x(1)^2);
        g(2) = -1/x(1)-1;
    case 'SkewNormal'
        m = par.snMu;
        Lambda = par.snSigma;
        alpha = par.snAlpha;
        g = -inv(Lambda)*(x-m)+mvnpdf(alpha'*(x-m)).*alpha./mvncdf(alpha'*(x-m));
    case 'MVN'
        mu = par.Mu;  Sigma= par.Sigma;
        g  = Sigma\(mu-x); 
   case 'MVNdiag'
        mu = par.Mu;  Sigma= par.Sigma;
        g  = diag(1./diag(Sigma))*(mu-x); 
end
end
