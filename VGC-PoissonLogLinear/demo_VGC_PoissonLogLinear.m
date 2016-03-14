% Poisson Log Linear Regression example for the paper
% "Variational Gaussian Copula Inference",
% Shaobo Han, Xuejun Liao, David. B. Dunson, and Lawrence Carin,
% The 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016)
% -----------------------------
% Code written by
% Shaobo Han, Duke University
% shaobohan@gmail.com
% 10/01/2015

% Methods
% 1 JAGS (Implemented via RJAGS)
% 2 Stochastic VGC with Bernstein Polynomials

% clear all; close all; clc;
load Treedata50.mat
% For more information about this dataset, see
% https://github.com/petrkeil/Statistics/blob/master/Lecture%203%20-%20poisson_regression/poisson_regression.Rmd

xm = x';
fix.n = length(y);
fix.Xmatrix = [ones(1, fix.n); xm; xm.^2]; % X is a 3 by N matrix
fix.y = y; fix.logfac_y =log(factorial(y));
fix.a0 = 1; fix.b0 = 1;
fix. P = 4; % Number of latent variables
figplot.nbins = 50;

trueModel = 'PoLog4'; opt.trueModel =trueModel;
VGmethod = 'Numeric';
inferModel2 = 'MVN'; % Multivariate Normal

%%  VGC Numeric
opt.adaptivePhi = 0;  opt.normalize = 1;
fix.c = 1;  % unnormalizing constant
opt.k = 10; % Degree/Maximum Degree of Bernstein Polynomials
opt.MaxIter = 50; % Number of SGD stages
opt.NumberZ = 1; % Average Gradients
opt.InnerIter = 500; % Number of iteration
opt.N_mc = 1; % Number of median average in sELBO stepsize search
opt.nsample = 1e5;

BPtype = 'BP';  % Bernstein Polynomials
% BPtype = 'exBP';  % Extended Bernstein Polynomials
PsiType = 'Normal';   opt.PsiPar(1) = 0;  opt.PsiPar(2) = 1;  % variance

opt.PoM.PsiType1 = 'Normal';
opt.PoM.PsiPar1(1) = 0;
opt.PoM.PsiPar1(2) = 1;  % variance
opt.PoM.PsiType2 = 'Exp';
opt.PoM.PsiPar2 = 1;

% PsiType = 'Exp';    opt.PsiPar = 1;
PhiType = 'Normal';  opt.PhiPar(1) = 0; opt.PhiPar(2) =1;  % variance
opt.Wthreshold = 1e4;

% Learning rate
opt.LearnRate.Mu = 1e-3; opt.LearnRate.C = 1e-4;
opt.LearnRate.W = 0.25.*1e-3; opt.LearnRate.dec = 0.95; % decreasing base learning rate

switch BPtype
    case 'BP'
        opt.D = opt.k;
    case 'exBP'
        opt.D = opt.k*(opt.k+1)/2; % # of basis functions
end

% Diagonal constraint on Upsilon
opt.diagUpsilon = 0;

%% Initialization
ini.Mu = opt.PhiPar(1).*ones(fix.P,1);
ini.C = sqrt(opt.PhiPar(2)).*0.1.*eye(fix.P);
% ini.w = randBPw(fix.P, opt.D, 1, 1);
ini.w = 1./opt.D.*ones(fix.P, opt.D);
% ini.w = randBPw(fix.P, opt.D, 1, 1);

% Median Outlier Removal
opt.OutlierTol = 100; % Threshold for online outlier detection
opt.WinSize = 20; % Size of the window
ini.WinSet = ini_WinSet(trueModel, PhiType, PsiType, BPtype, fix, opt, ini);

%% VGC MuCW-Numeric
opt.updateW = 1; opt.updateMuC = 1;
[ELBO3, par3] = vgcbp_MuCW(trueModel, inferModel2, PhiType, PsiType, BPtype, VGmethod, fix, ini, opt);

tic;
figplot.Y_svc3 = sampleGC(PhiType, PsiType, BPtype, opt, par3);
figplot.Y_svc3(sum(~isfinite(figplot.Y_svc3),2)~=0,:) = [];
toc;

PoissonLogVGC = figplot.Y_svc3';
save PoissonLogVGC.mat PoissonLogVGC
%% ------------------------------
load PoissonLogVGC.mat
load PoissonLogMCMC.mat

beta0 = beta0(:)'; beta1 = beta1(:)';
beta2 = beta2(:)'; tau = tau(:)';

X_mcmc = [beta0; beta1; beta2; tau];
X_vgc =PoissonLogVGC;
ncoutour = 6; figplot.ls = 1.5; figplot.nbins = 50;

figure
for i = 1:4
    for j = 1:4
        k = (i-1)*4+j;
        subplot(4,4,k)
        if i==j
            [f_x1,x_x1] = hist(X_mcmc(i,:), figplot.nbins);
            [f_x2,x_x2] = hist(X_vgc(i,:), figplot.nbins);
            plot(x_x1, f_x1./trapz(x_x1,f_x1), 'r--', 'linewidth', 2)
            hold on
            plot(x_x2, f_x2./trapz(x_x2,f_x2), 'b-', 'linewidth', 2)
            hold off
            legend('JAGS','VGC-BP')
        else
            tmp = hist1D2D([X_mcmc(i,:); X_mcmc(j,:)]', figplot.nbins);
            tmp2 = hist1D2D([X_vgc(i,:); X_vgc(j,:)]', figplot.nbins);
            contour(tmp.Seqx,tmp.Seqy, tmp.count./sum(sum(tmp.count)), ncoutour, 'r--' ,   'linewidth', 2*figplot.ls);
            hold on
            contour(tmp2.Seqx,tmp2.Seqy, tmp2.count./sum(sum(tmp2.count)), ncoutour, 'b' ,   'linewidth', 2*figplot.ls);
            hold off
            legend('JAGS', 'VGC-BP')
        end
    end
end
