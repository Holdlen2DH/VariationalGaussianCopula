function   [Delta_Mu, Delta_C, WinSet] = Grad_MuC(VGmethod, trueModel, inferModel, PhiType, PsiType, BPtype, fix, opt, par, WinSet)
NumberZ = opt.NumberZ;  P = fix.P;
PhiPar = opt.PhiPar; PsiPar = opt.PsiPar; k = opt.k;
Mu_old = par.Mu; C_old = par.C;  W_old = par.W;
MuMC = zeros(P, NumberZ);
CMC = zeros(P,P, NumberZ);
v_t = sqrt(diag(C_old*C_old'));

for jj = 1:NumberZ
    %%  Draw Sample w_t and Z_t
    Flag_outlier = 1;
    while Flag_outlier == 1
        % Generate a P by 1 multivariate Normal vector
        w_t = randn(P,1); Z_t = Mu_old + C_old*w_t;
        if opt.adaptivePhi == 1
            Z = diag(1./v_t)*(Z_t-Mu_old);
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
    if opt.adaptivePhi == 1
        sphi =  diag(1./v_t)*sphi_f(Z, PhiPar, PhiType);
        dphi = -diag(1./v_t)*(Z.*sphi);
    else
        sphi = sphi_f(Z_t, PhiPar, PhiType);
        dphi = dphi_f(Z_t, PhiPar, PhiType);
    end
    
    BPpdf_basiseval = BPpdf_basis(u, k, BPtype);
    tmpb = sum(BPpdf_basiseval.*W_old, 2);
    BPprime_basiseval = BPprime_basis(u, k, BPtype);
    tmpbprime = sum(BPprime_basiseval.*W_old, 2);
    rho1prime = sphi.*tmpbprime;
    
    spsi = spsi_f(tmpx, PsiPar, PsiType);
    dpsi = dpsi_f(tmpx,PsiPar, PsiType);
    g = derivatives(tmpx, trueModel, fix);
    
    hdz1 = diag(sphi./spsi)*tmpb;
    rho3prime = hdz1.*dpsi;
    hdz2 = (rho1prime.*sphi.*spsi...
        +tmpb.*dphi.*spsi...
        -tmpb.*sphi.*rho3prime)./(spsi.^2);
    lnh1dz =  hdz2./hdz1;
    ls_z = g.*hdz1+lnh1dz;
    switch VGmethod
        case 'Analytic'
            dmu = ls_z;
            dC = tril(ls_z*w_t')+ diag(1./diag(C_old));
        case 'Numeric'
            tmppar.Mu = Mu_old; tmppar.Sigma = C_old*C_old';
            g2 = logqpost_derivatives(Z_t, inferModel, tmppar);
            dmu = ls_z-g2;
            dC = tril(dmu*w_t');
    end
    MuMC(:, jj) = dmu;
    CMC(:, :, jj) = dC;
end
Delta_Mu = mean(MuMC, 2);
Delta_C= mean(CMC, 3);
end

