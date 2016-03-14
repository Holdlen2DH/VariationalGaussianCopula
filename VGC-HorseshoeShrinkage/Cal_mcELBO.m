function     sELBO = Cal_mcELBO(trueModel, inferModel, PhiType, PsiType, BPtype, fix, opt, par)
tmpMu=par.Mu;  tmpC = par.C;  tmpW=par.W;
N_mc = opt.N_mc; PhiPar = opt.PhiPar; PsiPar = opt.PsiPar; k = opt.k;
P = length(tmpMu); sELBOM = zeros(1, N_mc);
v_t = sqrt(diag(tmpC*tmpC'));

for i = 1: N_mc
    % Generate a P by 1 multivariate Normal vector
    w_t = randn(P,1);     Z_t = tmpMu + tmpC*w_t;
    if opt.adaptivePhi == 1
        Z = diag(1./v_t)*(Z_t-tmpMu);
        u = Phi_f(Z, PhiPar, PhiType);
        sphi = diag(1./v_t)*sphi_f(Z, PhiPar, PhiType);
    else
        u = Phi_f(Z_t, PhiPar, PhiType);
        sphi = sphi_f(Z_t, PhiPar, PhiType);
    end
    BPcdf_basiseval = BPcdf_basis(u, k, BPtype);
    tmpB = sum(BPcdf_basiseval.*tmpW, 2);
    tmpx = iPsi_f(tmpB, PsiPar, PsiType);
    logp = logmodel(tmpx, trueModel, fix);
    BPpdf_basiseval = BPpdf_basis(u, k, BPtype);
    tmpb = sum(BPpdf_basiseval.*tmpW, 2);
    spsi = spsi_f(tmpx, PsiPar, PsiType);
    hdz1 = diag(sphi./spsi)*tmpb;
    sELBOM(i) = logp + sum(log(hdz1))+ P/2*log(2*pi)+P/2 +  sum(log(diag(tmpC)));
end

% if opt.sELBO ==1
sELBO = median(sELBOM);
% % val = (P/2)*log(2*pi)+logdet(tmpC)+sELBO; % negative ELBO
% else
% sELBO+ P/2*log(2*pi)+P/2 +  sum(log(diag(C)));
% end

end
