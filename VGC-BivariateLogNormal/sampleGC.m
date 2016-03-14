function Y_svc = sampleGC(PhiType, PsiType, BPtype, opt, par)
PhiPar = opt.PhiPar; PsiPar = opt.PsiPar;
nsample = opt.nsample; k = opt.k;
Mu = par.Mu; C = par.C; W = par.W;
P = length(Mu);
% Sample from Gaussian Copula
optmu_s2 = Mu; optSig_s2 = C*C';
wtz2 = mvnrnd(optmu_s2, optSig_s2, nsample);
Y_svc = zeros(nsample,P);
if opt.adaptivePhi == 0
    for ns = 1:nsample
        tmpz = wtz2(ns,:)';
        tmpu = Phi_f(tmpz, PhiPar, PhiType); 
        BPcdf_basiseval = BPcdf_basis(tmpu, k, BPtype); 
        tmpB = sum(BPcdf_basiseval.*W, 2); 
        tmpx = iPsi_f(tmpB, PsiPar, PsiType); 
        Y_svc(ns,:) = tmpx;
    end
else    
    v_t = sqrt(diag(C*C'));
    for ns = 1:nsample
        tmpz0 = wtz2(ns,:)';
        tmpz = diag(1./v_t)*(tmpz0-Mu);
        tmpu = Phi_f(tmpz, PhiPar, PhiType); 
        BPcdf_basiseval = BPcdf_basis(tmpu, k, BPtype); 
        tmpB = sum(BPcdf_basiseval.*W, 2); 
        tmpx = iPsi_f(tmpB, PsiPar, PsiType); 
        Y_svc(ns,:) = tmpx;
    end
end
end