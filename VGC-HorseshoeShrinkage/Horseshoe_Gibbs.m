function Y_mc=Horseshoe_Gibbs(x, nsample)
tic; 
Y_mc = zeros(nsample,2); 
for iter = 1:nsample
    if iter==1; ga_mc = 1; end
    tau_mc = 1/gamrnd(1,1/(x^2/2+ga_mc)); 
    ga_mc = gamrnd(1,1/(1/tau_mc+1)); 
    Y_mc(iter,:) = [tau_mc; ga_mc]; 
end
t2=toc; 
display(['[Gibbs Sampler] Num of Iters: ' num2str(nsample),' Elapsed Time: ', num2str(t2) ' sec']);
tempC=corrcoef(log(Y_mc)); 
display(['[Gibbs Sampler] Correlation Coefs:  C11 = ', num2str(tempC(1,1)), ...
    ', C21 = ', num2str(tempC(2,1)), ', C22 =  ', num2str(tempC(2,2))])
end