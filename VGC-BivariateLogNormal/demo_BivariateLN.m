% Bivariate Log-Normal example for the paper 
% "Variational Gaussian Copula Inference", 
% Shaobo Han, Xuejun Liao, David. B. Dunson, and Lawrence Carin, 
% The 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016)
% -----------------------------
% Code written by 
% Shaobo Han, Duke University
% shaobohan@gmail.com

% 09/23/2015

% Method
% 1 Ground Truth
% 2  VGC-LN with Analytic VG updates
% 3  VGC-LN with Numeric VG updates
% 4  VGC-BP with Numeric VG updates
% 5  VGC-BP with Analytic VG updates

% We use MVLOGNRAND written by Stephen Lienhard 
% as Random Samples Generator for Bivaraite Lognormal Distributions 
% http://www.mathworks.com/matlabcentral/fileexchange/6426-multivariate-lognormal-simulation-with-correlation/content/MvLogNRand.m

clear all; close all; clc;
% addpath(genpath(pwd))     % Add all sub-directories to the Matlab path
fix. P = 2; % Number of latent variables

%  True SN Posterior
trueModel = 'LogNormal2';
fix.LNsigmaY1 = 0.5; fix.LNsigmaY2 = 0.5;
fix.LNmuY1 = 0.1;  fix.LNmuY2 = 0.1;  
fix.LNrho = -0.4; % fix.LNrho = 0.4;
 
 fix.trueMu = [fix.LNmuY1;  fix.LNmuY2];
 fix.trueSigma=[fix.LNsigmaY1^2,fix.LNrho*fix.LNsigmaY1*fix.LNsigmaY2;...
     fix.LNrho*fix.LNsigmaY1*fix.LNsigmaY2,fix.LNsigmaY2^2];
fix.trueUpsilon = corrcov( fix.trueSigma);
opt.adaptivePhi = 0;  opt.normalize = 0; 
fix.c = 2;  % unnormalizing constant 
%%
opt.k = 5; % Degree/Maximum Degree of Bernstein Polynomials
opt.MaxIter = 100; % Number of SGD stages
opt.NumberZ = 1; % Average Gradients
opt.InnerIter = 250; % Number of iteration

stageMax = 100; iterMax = 25;
fix.rho = 0.005; fix.dec = 0.95;

opt.N_mc = 1; % Number of median average in sELBO stepsize search
opt.nsample = 5e5;
figplot.nbins = 50;
BPtype = 'BP';  % Bernstein Polynomials
% BPtype = 'exBP';  % Extended Bernstein Polynomials

% PsiType = 'Normal';   opt.PsiPar(1) = 0;  opt.PsiPar(2) = 1;  % variance
PsiType = 'Exp';    opt.PsiPar = 1;
PhiType = 'Normal';  opt.PhiPar(1) = 0; opt.PhiPar(2) =1;  % variance
VGmethod = 'Numeric'; 
opt.Wthreshold = 1e12; 

% Learning rate
opt.LearnRate.Mu = 0.005; opt.LearnRate.C = 0.005;
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
ini.C = sqrt(opt.PhiPar(2)).*eye(fix.P);
% ini.w = randBPw(fix.P, opt.D, 1, 1);
ini.w = 1./opt.D.*ones(fix.P, opt.D);
% ini.w = randBPw(fix.P, opt.D, 1, 1);
 
% Median Outlier Removal
opt.OutlierTol = 10; % Threshold for online outlier detection
opt.WinSize = 20; % Size of the window
ini.WinSet = ini_WinSet(trueModel, PhiType, PsiType, BPtype, fix, opt, ini);

%%  1 VGC LN Analytic
fix.scale = 1; % The scale of covariance
inferModel2 = 'MVN'; % Multivariate Normal
VG = svgln(trueModel, inferModel2, stageMax, iterMax, 'Analytic', fix);
figplot.Y_vg = sampleEVG(opt, VG); 
figplot.Y_vg(sum(~isfinite(figplot.Y_vg),2)~=0,:) = [];
VGsamples = hist1D2D(figplot.Y_vg, figplot.nbins);

%%  VG Numeric
VG2 = svgln(trueModel, inferModel2, stageMax, iterMax, 'Numeric', fix);
figplot.Y_vg2 = sampleEVG(opt, VG2); 
figplot.Y_vg2(sum(~isfinite(figplot.Y_vg2),2)~=0,:) = [];
VGsamples2 = hist1D2D(figplot.Y_vg2, figplot.nbins);

%% VGC MuCW-Numeric
opt.updateW = 1; opt.updateMuC = 1; 
[ELBO3, par3] = vgcbp_MuCW(trueModel, inferModel2, PhiType, PsiType, BPtype, 'Numeric', fix, ini, opt); 
% ini.Mu = par.Mu; ini.C = par.C; ini.WinSet = par.WinSet; 

figplot.Y_svc3 = sampleGC(PhiType, PsiType, BPtype, opt, par3);
figplot.Y_svc3(sum(~isfinite(figplot.Y_svc3),2)~=0,:) = [];
VGCBPall3 = hist1D2D(figplot.Y_svc3, figplot.nbins);

%% VGC MuCW-Analytic
opt.updateW = 1; opt.updateMuC = 1; 
[ELBO5, par5] = vgcbp_MuCW(trueModel, inferModel2, PhiType, PsiType, BPtype, 'Analytic', fix, ini, opt); 
figplot.Y_svc5 = sampleGC(PhiType, PsiType, BPtype, opt, par5);
figplot.Y_svc5(sum(~isfinite(figplot.Y_svc5),2)~=0,:) = [];
VGCBPall5 = hist1D2D(figplot.Y_svc5, figplot.nbins);

%%  Plot
figplot.nb = 1000;  
contourvector = [0.05, 0.1, 0.25,0.5, 0.75, 0.9]; 
figplot.Seqx = linspace(0.001, 5, figplot.nb)';
figplot.Seqy = linspace(0.001, 5, figplot.nb+1)';
[figplot.X1,figplot.X2] = meshgrid(figplot.Seqx,figplot.Seqy);
figplot.log_pM = modelLN(figplot.X1, figplot.X2, fix);
plevs1=contourvector.*max(exp(figplot.log_pM (:)));

figure (1)
figplot.ls = 2; lfs =26; 
set(gca, 'fontsize', 16)
contour(figplot.Seqx,figplot.Seqy, exp(figplot.log_pM), plevs1,  'linewidth', 2*figplot.ls);
xlabel('x_1'); ylabel('x_2');
axis([0,3.5,0,3.5])
AX=legend('Ground Truth');
LEG = findobj(AX,'type','text');
set(LEG,'fontsize',lfs)

figure (2)
F25 =VGsamples.count./sum(sum(VGsamples.count)); 
plevs3=contourvector.*max(F25(:));
set(gca, 'fontsize', 16)
contour(VGsamples.Seqx,VGsamples.Seqy, F25, plevs3,  'linewidth', 2*figplot.ls);
xlabel('x_1'); ylabel('x_2');
axis([0,3.5,0,3.5])
AX= legend('VGC-LN1');
LEG = findobj(AX,'type','text');
set(LEG,'fontsize',lfs)
 
figure (3)
F26 =VGsamples2.count./sum(sum(VGsamples2.count)); 
plevs4=contourvector.*max(F26(:));
set(gca, 'fontsize', 16)
contour(VGsamples2.Seqx,VGsamples2.Seqy, F26, plevs4,  'linewidth', 2*figplot.ls);
xlabel('x_1'); ylabel('x_2');
axis([0,3.5,0,3.5])
AX= legend('VGC-LN2');
LEG = findobj(AX,'type','text');
set(LEG,'fontsize',lfs)

figure (4)
set(gca, 'fontsize', 16)
F38 =VGCBPall5.count./sum(sum(VGCBPall5.count)); 
plevs36=contourvector.*max(F38(:));
contour(VGCBPall5.Seqx,VGCBPall5.Seqy, F38, plevs36,  'linewidth', 2*figplot.ls);
xlabel('x_1'); ylabel('x_2');
axis([0,3.5,0,3.5])
AX=legend('VGC-BP1');
LEG = findobj(AX,'type','text');
set(LEG,'fontsize',lfs)

figure (5)
set(gca, 'fontsize', 16)
F28 =VGCBPall3.count./sum(sum(VGCBPall3.count)); 
plevs6=contourvector.*max(F28(:));
contour(VGCBPall3.Seqx,VGCBPall3.Seqy, F28, plevs6,  'linewidth', 2*figplot.ls);
xlabel('x_1'); ylabel('x_2');
axis([0,3.5,0,3.5])
AX=legend('VGC-BP2');
LEG = findobj(AX,'type','text');
set(LEG,'fontsize',lfs)

figure (6)
figplot.fs = 16; figplot.ms =3;
set(gca, 'fontsize', figplot.fs)
plot(median(VG.RMSECTC,2), '--r*', 'linewidth', figplot.ls, 'Markersize', 3+figplot.ms);
hold on 
plot(median(VG2.RMSECTC,2), '-bo', 'linewidth', figplot.ls, 'Markersize', 3+figplot.ms);
hold off
grid on
axis([0,100, 0,1])
AX=legend('VGC-LN1', 'VGC-LN2');
LEG = findobj(AX,'type','text');
set(LEG,'fontsize',lfs)
xlabel('Iterations')
ylabel('\rho')

figure (7)
figplot.fs = 16; figplot.ms = 3;
set(gca, 'fontsize', figplot.fs)
plot([par5.RMSE(1,1), median(par5.RMSE,1)], '--r*', 'linewidth', figplot.ls, 'Markersize', 3+figplot.ms);
hold on 
plot([par3.RMSE(1,1),median(par3.RMSE,1)], '-bo', 'linewidth', figplot.ls, 'Markersize', 3+figplot.ms);
hold off
grid on
axis([0,100, 0,1])
AX=legend('VGC-BP1', 'VGC-BP2');
LEG = findobj(AX,'type','text');
set(LEG,'fontsize',lfs)
xlabel('Iterations')
ylabel('\rho')