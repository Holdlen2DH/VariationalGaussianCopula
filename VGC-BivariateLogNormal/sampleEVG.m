function Y_svc = sampleEVG(opt, VG)
nsample = opt.nsample; 
Mu = VG.mu; C = VG.C; 
% wtz2 = mvnrnd(Mu', C*C', nsample);
% Y_svc = exp(wtz2); 
Sigma = sqrt(diag(C*C')); Simulations = nsample; 
CorrMat = corrcov(C*C'); 
Y_svc = MvLogNRand( Mu , Sigma , Simulations , CorrMat ); 
end