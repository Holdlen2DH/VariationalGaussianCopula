function     val = Cal_mcELBO(trueModel, inferModel, PhiType, PsiType, BPtype, fix, opt, par)
tmpMu=par.Mu;  tmpC = par.C;  tmpW=par.W;
N_mc = opt.N_mc; PhiPar = opt.PhiPar; PsiPar = opt.PsiPar; k = opt.k;
P = length(tmpMu); sELBOM = zeros(1, N_mc);
for i = 1: N_mc
    % Generate a P by 1 multivariate Normal vector
    w_t = randn(P,1); Z_t = tmpMu + tmpC*w_t;
    u = Phi_f(Z_t, PhiPar, PhiType);
    BPcdf_basiseval = BPcdf_basis(u, k, BPtype);
    tmpB = sum(BPcdf_basiseval.*tmpW, 2); 
    tmpx = iPsi_f(tmpB, PsiPar, PsiType); 
    logp = logmodel(tmpx, trueModel, fix);
    BPpdf_basiseval = BPpdf_basis(u, k, BPtype); 
    tmpb = sum(BPpdf_basiseval.*tmpW, 2); 
    sphi = sphi_f(Z_t, PhiPar, PhiType);
    spsi = spsi_f(tmpx, PsiPar, PsiType);
    hdz1 = diag(sphi./spsi)*tmpb;
    logq = logqpost(Z_t, inferModel, par); 
    sELBOM(i) = logp + sum(log(hdz1))-logq;
end
sELBO = median(sELBOM);
val = sELBO; 
end