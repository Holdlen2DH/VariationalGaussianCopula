function   [Delta_w, WinSet] = Grad_W(trueModel, PhiType, PsiType, BPtype, fix, opt, par, WinSet)
NumberZ = opt.NumberZ;  P = fix.P; D = opt.D;
PhiPar = opt.PhiPar; PsiPar = opt.PsiPar; k = opt.k;
Mu_old = par.Mu; C_old = par.C; W_old = par.W;
WMC = zeros(P, D, NumberZ);
vt = sqrt(diag(C_old*C_old')); 
for jj = 1:NumberZ
    %%  Draw Sample w_t and Z_t
    Flag_outlier = 1;
    while Flag_outlier == 1
        % Generate a P by 1 multivariate Normal vector
        w_t = randn(P,1); Z_t = Mu_old + C_old*w_t;
        if opt.adaptivePhi == 1
        Z = diag(1./vt)*(Z_t-Mu_old); 
        u = Phi_f(Z, PhiPar, PhiType);
        else
        u = Phi_f(Z_t, PhiPar, PhiType);
        end
        BPcdf_basiseval = BPcdf_basis(u, k, BPtype); 
        tmpB = sum(BPcdf_basiseval.*W_old, 2); 
        tmpx = iPsi_f(tmpB, PsiPar, PsiType); 
        logp = logmodel(tmpx, trueModel, fix);
        % Online outlier detection
        score = mzscore_test(logp, WinSet);
        if abs(score) < opt.OutlierTol
            Flag_outlier = 0; WinSet_old = WinSet;
            WinSet = [WinSet_old, logp];
            WinSet(1) = [];
        end
    end
    BPpdf_basiseval = BPpdf_basis(u, k, BPtype); 
    tmpb = sum(BPpdf_basiseval.*W_old, 2); 
    spsi = spsi_f(tmpx, PsiPar, PsiType);
    dpsi = dpsi_f(tmpx,PsiPar, PsiType);
    g = derivatives(tmpx, trueModel, fix);
    hdw = diag(1./spsi)*BPcdf_basiseval;
    lnh1dw = diag(1./tmpb)*BPpdf_basiseval - diag(dpsi./spsi.^2)*BPcdf_basiseval;
   WMC(:,:,jj)  = diag(g)*hdw+lnh1dw;
end
Delta_w = median(WMC, 3);
end