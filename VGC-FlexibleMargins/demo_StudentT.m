% Flexible Margins example for the paper
% "Variational Gaussian Copula Inference",
% Shaobo Han, Xuejun Liao, David. B. Dunson, and Lawrence Carin,
% The 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016)
% -----------------------------
% Code written by
% Shaobo Han, Duke University
% shaobohan@gmail.com
% 09/22/2015

% Examples:
%  Student's t Distribution

clear all; close all; clc;
% addpath(genpath(pwd))     % Add all sub-directories to the Matlab path

fix. P = 1; % Number of latent variables
trueModel = 'Student';
fix.nu = 1;
inferModel ='MVNdiag';
fix.c = 1;  % unnormalizing constant
%%
opt.k =15; % Degree/Maximum Degree of Bernstein Polynomials
opt.MaxIter = 20; % Number of SGD stages
opt.NumberZ = 1; % Average Gradients
opt.InnerIter = 500; % Number of iteration
opt.N_mc = 1; % Number of median average in sELBO stepsize search

BPtype = 'BP';  % Bernstein Polynomials
% BPtype = 'exBP';  % Extended Bernstein Polynomials

PsiType = 'Normal';   opt.PsiPar(1) = 0;  opt.PsiPar(2) = 1;  % variance
% PsiType = 'Exp';    opt.PsiPar = 0.5;
PhiType = 'Normal';  opt.PhiPar(1) = 0; opt.PhiPar(2) =1;  % variance

% Learning rate
opt.LearnRate.Mu = 0.001; opt.LearnRate.C =1e-3;
opt.LearnRate.W = 0.5*1e-3; opt.LearnRate.dec = 0.95; % decreasing base learning rate

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
ini.C = opt.PhiPar(2).*eye(fix.P);
% ini.w = randBPw(fix.P, opt.D, 1, 1);
ini.w = 1./opt.D.*ones(fix.P, opt.D);

% Median Outlier Removal
opt.OutlierTol = 10; % Threshold for online outlier detection
opt.WinSize = 20; % Size of the window

opt.Wthreshold = 1e12;  opt.normalize = 0;
ini.WinSet = ini_WinSet(trueModel, PhiType, PsiType, BPtype, fix, opt, ini);
%% VGC-(Uniform)Adaptive Algorithm
[ELBO, par]   = vgcbp_w(trueModel, inferModel, PhiType, PsiType, BPtype,fix, ini, opt);
opt.nsample = 5e5;
Y_svc = sampleGC(PhiType, PsiType, BPtype, opt, par);
Y_svc(sum(~isfinite(Y_svc),2)~=0,:) = [];
figplot.nbins = 50;
[figplot.f_x1,figplot.x_x1] = hist(Y_svc, figplot.nbins);

%% VGC-(Non-Uniform) Algorithm
opt.adaptivePhi = 0; VGmethod = 'Numeric';
[ELBO2, par2]  = vgcbp_MuCW(trueModel, inferModel, PhiType, PsiType, BPtype, VGmethod, fix, ini, opt);
Y_svc2 = sampleGC(PhiType, PsiType, BPtype, opt, par2);
Y_svc2(sum(~isfinite(Y_svc2),2)~=0,:) = [];
[figplot.f_x2,figplot.x_x2] = hist(Y_svc2, figplot.nbins);

%%  Plot
figplot.nb = 100;  figplot.Seqx = linspace(-8, 8, figplot.nb)';
figplot.log_pM =  logmodel(figplot.Seqx, trueModel, fix);
trueModel2 = 'SkewNormal';
fix2.snMu = 0;  fix2.snSigma = 1; fix2.snAlpha =0; fix2.c = 1;
figplot.log_pM2 =  logmodel(figplot.Seqx, trueModel2, fix2);

figure
figplot.ls =1; figplot.ms =1;
set(gca,'fontsize', 16)
plot(figplot.Seqx, exp(figplot.log_pM)./trapz(figplot.Seqx, exp(figplot.log_pM)), '-', 'linewidth', 3.5*figplot.ls, 'Markersize', 4*figplot.ms, 'Color', 'k');
hold on
plot(figplot.x_x1,figplot.f_x1/trapz(figplot.x_x1,figplot.f_x1),  '--', 'linewidth', 2.5*figplot.ls, 'Markersize', 4*figplot.ms, 'Color', 'g');
hold on
plot(figplot.x_x2,figplot.f_x2/trapz(figplot.x_x2,figplot.f_x2), '-o', 'linewidth', 2*figplot.ls, 'Markersize', 4*figplot.ms, 'Color', 'b');
hold on
plot(figplot.Seqx, exp(figplot.log_pM2)./trapz(figplot.Seqx, exp(figplot.log_pM2)), '-.', 'linewidth', 4*figplot.ls, 'Markersize', 4*figplot.ms, 'Color', 'r');
hold off
legend('Ground Truth', 'VIT-BP(k=15)', 'VGC-BP(k=15)' , 'Predefined \Psi')

