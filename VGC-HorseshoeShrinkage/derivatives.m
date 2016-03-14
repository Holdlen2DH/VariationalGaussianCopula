% Derivatives of Log likelihood and Log prior
% Model dependent term
% Shaobo Han
% 08/10/2015

function g = derivatives(x, model, fix)
switch model
    case 'Horseshoe'
        y = fix.y;
        g = zeros(2,1);
        g(1) = -2/x(1)+y^2/2/(x(1)^2)+x(2)/(x(1)^2);
        g(2) = -1/x(1)-1; 
    case 'SkewNormal'
        m = fix.snMu;
        Lambda = fix.snSigma;
        alpha = fix.snAlpha;
        tmp = alpha'*(x-m);
        g = -(Lambda\(x-m))+normpdf(tmp).*alpha./normcdf(tmp);
    case 'LogNormal2'
        sigmaY1 = fix.LNsigmaY1; sigmaY2 = fix.LNsigmaY2;
        muY1 = fix.LNmuY1; muY2 = fix.LNmuY2;  LNrho = fix.LNrho;
        alpha1 = (log(x(1))-muY1)./sigmaY1;
        alpha2= (log(x(2))-muY2)./sigmaY2;
        g = zeros(2,1);
        g(1)=-1./x(1)-(alpha1-LNrho*alpha2)/(1-LNrho.^2)./x(1)./sigmaY1; 
        g(2)=-1./x(2)-(alpha2-LNrho*alpha1)/(1-LNrho.^2)./x(2)./sigmaY2; 
    case 'Gamma'
        alpha = fix.alpha;
        beta = fix.beta;
        g = (alpha-1)./x-beta;
    case 'Student'
        nu = fix.nu;
        g = -(nu+1)*x/(nu+x^2);
    case 'Beta';
        alpha = fix.alpha;  beta = fix.beta;
        g = (alpha-1)/x-(beta-1)/(1-x);
end
end
