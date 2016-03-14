function [Y_vb,ELBO]=Horseshoe_Mfvb(x,nsample, iterMax,tol)
tic;
c0 = -0.5*log(2*pi)-2*log(gamma(0.5));
iter = 1;  ga_exp = 1;
tau_a = 1; ga_a = 1; Flag=1;
ELBO = zeros(1, iterMax);
while (iter<=iterMax) && (Flag==1)
    tau_b = x^2/2+ga_exp;
    invtau_exp = tau_a/tau_b;
    ga_b = invtau_exp+1;
    ga_exp = ga_a/ga_b;
    H1 = tau_a+log(tau_b)+log(gamma(tau_a))-(1+tau_a)*psi(tau_a);
    H2 = ga_a-log(ga_b)+log(gamma(ga_a))+(1-ga_a)*psi(ga_a);
    ELBO(iter) = c0-2*(log(tau_b)-psi(tau_a))-x^2/2*tau_a/tau_b-1+H1+H2;
    if iter>1
        if norm(ELBO(iter)-ELBO(iter-1),'fro')/norm(ELBO(iter),'fro')<tol
            Flag = 0;
        end
        if ELBO(iter)<ELBO(iter-1)
            display('Error: ELBO is decreasing!!')
        end
    end
    iter = iter+1;
end
t2=toc;
ELBO = ELBO(1:(iter-1));
tau_mfvb = 1./gamrnd(tau_a, 1/tau_b, nsample,1);
ga_mfvb = gamrnd(ga_a, 1/ga_b, nsample,1);
Y_vb = [tau_mfvb, ga_mfvb];
display(['[Mean Field VB] ELBO = ', num2str(ELBO(iter-1)),  ', Num of Iters: ' num2str(iter-1),' Elapsed Time: ', num2str(t2) ' sec']);
tempC=corrcoef(log(Y_vb)); 
display(['[Mean Field VB] Correlation Coefs:  C11 = ', num2str(tempC(1,1)), ...
    ', C21 = ', num2str(tempC(2,1)), ', C22 =  ', num2str(tempC(2,2))])
end