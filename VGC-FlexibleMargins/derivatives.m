% Derivatives of Log likelihood and Log prior
% Model dependent term
% Shaobo Han
% 08/10/2015

function g = derivatives(x, model, fix)
switch model
    case 'SkewNormal'
        m = fix.snMu; 
        Lambda = fix.snSigma; 
        alpha = fix.snAlpha; 
        g = -(Lambda\(x-m))+mvnpdf(alpha'*(x-m)).*alpha./mvncdf(alpha'*(x-m)); 
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
