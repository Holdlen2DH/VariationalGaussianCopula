% Log likelihood and Log prior
% Model dependent term
% Shaobo Han
% 08/10/2015

function p =  logmodel(x, model, fix)
switch model
    case 'PoLog4'
        y = fix.y; Xmatrix = fix.Xmatrix;   logfac_y=fix.logfac_y;
        mu = exp(Xmatrix'*x(1:3));
        tmp0 = sum(y.*log(mu)-mu-logfac_y);
        tmp_tau = x(4); tmp_sqrttau = sqrt(tmp_tau);
        tmp1 = log(normpdf(x(1), 0, tmp_sqrttau));
        tmp2 = log(normpdf(x(2), 0, tmp_sqrttau));
        tmp3 = log(normpdf(x(3), 0, tmp_sqrttau));
        tmp4 = log(gampdf(x(4), fix.a0, 1./fix.b0));
        p = tmp0+tmp1+tmp2+tmp3+tmp4;
    case 'Horseshoe'
        y = fix.y;
        c0 = -0.5*log(2*pi)-2*log(gamma(0.5));
        p = c0-2*log(x(1))-y.^2./2./x(1)-x(2)./x(1)-x(2);
    case 'LogNormal2'
        sigmaY1 = fix.LNsigmaY1; sigmaY2 = fix.LNsigmaY2;
        muY1 = fix.LNmuY1; muY2 = fix.LNmuY2;  LNrho = fix.LNrho;
        c0 = -log(2*pi)-log(sigmaY1*sigmaY2*sqrt(1-LNrho^2));
        alpha1 = (log(x(1))-muY1)/sigmaY1;
        alpha2= (log(x(2))-muY2)/sigmaY2;
        q = 1/(1-LNrho^2)*(alpha1^2-2*LNrho*alpha1*alpha2+alpha2^2);
        p = c0-log(x(1))-log(x(2))-q/2;
    case 'LogNormal'
        sigmaY = fix.LNsigmaY;
        muY = fix.LNmuY;
        c0 = -0.5*log(2*pi)-log(sigmaY);
        alpha = (log(x)-muY)./sigmaY;
        p = c0-log(x)- alpha.^2/2;
    case 'SkewNormal'
        m = fix.snMu;
        Lambda = fix.snSigma;
        alpha = fix.snAlpha;
        c = fix.c;
        p = log(2)+log(mvnpdf(x,m,Lambda))+log(normcdf(alpha'*(x-m)))+log(c);
    case 'Gamma'
        alpha = fix.alpha;
        beta = fix.beta;
        c = fix.c;
        p = log(gampdf(x, alpha, 1/beta))+log(c);
    case 'Student'
        nu = fix.nu;         c = fix.c;
        p = log(tpdf(x,nu))+log(c);
    case 'Exp'
        lambda = fix.lambda;         c = fix.c;
        p = log(exppdf(x,1/lambda))+log(c);
    case 'Beta';
        alpha = fix.alpha;  beta = fix.beta;  c = fix.c;
        p = log(betapdf(x, alpha, beta))+log(c);
end
end
