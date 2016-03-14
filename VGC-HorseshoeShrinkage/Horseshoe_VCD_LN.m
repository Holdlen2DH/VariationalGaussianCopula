function [Y_vc,ELBO,postoptimized]=Horseshoe_VCD_LN(post, options,nsample, fix)
%  Log-normal margins (exponential transform)
tic;
if fix.Flagdiag ==1
    post.c(2,1) = 0;
    [xmin,f,~,output] =  minFunc(@vcval_Horseshoe_diag,post2vec(post),options,fix);
else
    [xmin,f,~,output] =  minFunc(@vcval_Horseshoe,post2vec(post),options,fix);
end
postoptimized = vec2post(xmin,fix.post.dim);
t2 = toc;
optmu = postoptimized.m;
optSig = postoptimized.c*postoptimized.c';
wtz = mvnrnd(optmu, optSig, nsample);
tau_vc = exp(wtz(:,1));
gamma_vc = exp(wtz(:,2));
Y_vc = [tau_vc,  gamma_vc];
ELBO=-output.trace.fval;
tempC=corrcoef(log(Y_vc));

if fix.Flagdiag ==1
    display(['[VC(D): LogNormal+diag] ELBO = ', num2str(-f),  ', Num of Iters: ' num2str(output.iterations),' Elapsed Time: ', num2str(t2) ' sec']);
    display(['[VC(D): LogNormal+diag] Correlation Coefs:  C11 = ', num2str(tempC(1,1)), ...
        ', C21 = ', num2str(tempC(2,1)), ', C22 =  ', num2str(tempC(2,2))])
else
    display(['[VC(D): LogNormal+full] ELBO = ', num2str(-f),  ', Num of Iters: ' num2str(output.iterations),' Elapsed Time: ', num2str(t2) ' sec']);
    display(['[VC(D): LogNormal+full] Correlation Coefs:  C11 = ', num2str(tempC(1,1)), ...
        ', C21 = ', num2str(tempC(2,1)), ', C22 =  ', num2str(tempC(2,2))])
end
end