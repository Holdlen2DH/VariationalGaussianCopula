function WinSet = ini_WinSet(trueModel, PhiType, PsiType, BPtype, fix, opt, ini)
k = opt.k;  PhiPar = opt.PhiPar; PsiPar = opt.PsiPar;
P = fix.P;  Mu_old = ini.Mu; C_old = ini.C; W_old = ini.w;
WinSet = zeros(1, opt.WinSize);
for i = 1:opt.WinSize
    w_t0 = randn(P,1); Z_t0 = Mu_old + C_old*w_t0;
    u0 = Phi_f(Z_t0, PhiPar, PhiType);
    BPcdf_basiseval0 = BPcdf_basis(u0, k, BPtype); 
    tmpB0 = sum(BPcdf_basiseval0.*W_old, 2); 
    tmpx0 = iPsi_f(tmpB0, PsiPar, PsiType);
    WinSet(i) = logmodel(tmpx0, trueModel, fix);
end
end