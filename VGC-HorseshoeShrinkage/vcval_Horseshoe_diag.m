function [val, grad] = vcval_NIGG_diag(opt,fix)
% Evaluate ELBO and Gradient  
optl = length(opt);                             % get size of opt params
% break if opt vals not real
if ~isreal(opt), val = Inf; grad = Inf(optl,1); return; end
post  =  vec2post(opt,fix.post.dim); % unpack opt. params.
% unpack fix. parameters 
c1  =  fix.c1; 
x  =  fix.x; 

%% ELBO
% ENTROPY TERM
Eh=log(det(post.c));
m1 = post.m(1); m2 = post.m(2); 
C11 = post.c(1,1); C21 = 0; C22 = post.c(2,2); 
% break if diag(C) negative
if ~isreal(Eh), val = Inf; grad = Inf(optl,1); return; end
ell1 = exp(-m1+(C11^2)/2); 
tempC = C11^2-2*C11*C21+C21^2+C22^2; 
ell2 = exp(m2-m1+0.5.*tempC); 
ell3 = exp(m2+0.5*(C21^2+C22^2)); 
val1  =  c1-m1+m2-x^2/2*ell1-ell2-ell3; 
val2  =  Eh; 
val  =  -(val1+val2); 
if ~isreal(val), val = Inf; grad = Inf(optl,1); return; end
%% Gradients 
grad_m1  =  -1+x^2/2*ell1+ell2; 
grad_m2  =  1-ell2-ell3; 
tempgrad_c11  =  -C11*(x^2)/2*ell1-(C11-C21)*ell2+1/C11; 
tempgrad_c21  =  0;
tempgrad_c22  =  -C22*ell2-C22*ell3+1/C22; 
gradc  =  [tempgrad_c11, 0; tempgrad_c21, tempgrad_c22]; 
gradm  =  [grad_m1; grad_m2]; 
grad  =  [-gradm; -gradc(:)];
end