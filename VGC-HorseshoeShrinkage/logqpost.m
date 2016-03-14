% Log q posteiror 
% Model dependent term
% Shaobo Han
% 09/11/2015

function p =  logqpost(x, model, par)
switch model
    case 'Horseshoe'
        y = par.y; 
        c0 = -0.5*log(2*pi)-2*log(gamma(0.5));
        p = c0-2*log(x(1))-y.^2./2./x(1)-x(2)./x(1)-x(2);
    case 'SkewNormal'        
        m = par.snMu; 
        Lambda = par.snSigma; 
        alpha = par.snAlpha; 
        p = log(2)+log(mvnpdf(x,m,Lambda))+log(normcdf(alpha'*(x-m))) ;
    case 'MVN'
        mu = par.Mu;  C= par.C; Sigma=C*C';
        d = length(mu); 
        p = -d*log(2*pi)/2-0.5*logdet(Sigma)-0.5*(x-mu)'*(Sigma\(x-mu)); 
     case 'MVNdiag'
        mu = par.Mu;  C= par.C; Sigma=C*C';
        p = sum(log(normpdf(x, mu, sqrt(diag(Sigma))))); 
end
end
