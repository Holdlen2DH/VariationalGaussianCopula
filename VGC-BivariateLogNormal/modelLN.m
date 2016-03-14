function p = modelLN(pi_1V, pi_2V, fix)
[m,n] = size(pi_1V); 
sigmaY1 = fix.LNsigmaY1; sigmaY2 = fix.LNsigmaY2; 
muY1 = fix.LNmuY1; muY2 = fix.LNmuY2;  LNrho = fix.LNrho; 
x1 = reshape(pi_1V, m*n,1); 
x2 = reshape(pi_2V, m*n,1); 
c0 = -log(2*pi)-log(sigmaY1*sigmaY2*sqrt(1-LNrho^2)); 
alpha1 = (log(x1)-muY1)./sigmaY1; 
alpha2 = (log(x2)-muY2)./sigmaY2; 
q = 1/(1-LNrho^2).*(alpha1.^2-2.*LNrho.*alpha1.*alpha2+alpha2.^2); 
logpV = c0-log(x1)-log(x2)-q/2; 
 p = reshape(logpV, m,n); 
end
