% Log likelihood and Log prior
% Model dependent term
% Shaobo Han
% 08/10/2015

function p =  logmodel(x, model, fix)
switch model
    case 'SkewNormal'
        m = fix.snMu;
        Lambda = fix.snSigma;
        alpha = fix.snAlpha;
        c = fix.c;
        p = log(2)+log(mvnpdf(x,m,Lambda))+log(normcdf(alpha'*(x-m)))+log(c) ;
    case 'Gamma'
        alpha = fix.alpha;
        beta = fix.beta;
        c = fix.c;
        % p = alpha*log(beta)-log(gamma(alpha))+(alpha-1)*log(x)-beta*x+log(c) ;
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
