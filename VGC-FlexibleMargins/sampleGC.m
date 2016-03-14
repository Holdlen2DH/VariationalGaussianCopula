function Y_svc = sampleGC(PhiType, PsiType, BPtype, opt, par)
PhiPar = opt.PhiPar; PsiPar = opt.PsiPar;
nsample = opt.nsample; k = opt.k;
Mu = par.Mu; C = par.C; W = par.W;
% Sample from Gaussian Copula
optmu_s2 = Mu; optSig_s2 = C*C';
wtz2 = mvnrnd(optmu_s2, optSig_s2, nsample);
tau_svc2 = zeros(nsample,1);
for ns = 1:nsample
    tmpz = wtz2(ns,:)';
    tmpu = Phi_f(tmpz, PhiPar, PhiType); 
    BPcdf_basiseval = BPcdf_basis(tmpu, k, BPtype); 
    tmpB = sum(BPcdf_basiseval.*W, 2); 
    tmpx = iPsi_f(tmpB, PsiPar, PsiType); 
    tau_svc2(ns,1) = tmpx(1);
end
Y_svc = tau_svc2;
end