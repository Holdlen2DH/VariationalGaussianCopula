% Horseshoe shrinakge example for the paper
% "Variational Gaussian Copula Inference",
% Shaobo Han, Xuejun Liao, David. B. Dunson, and Lawrence Carin,
% The 19th International Conference on Artificial Intelligence and Statistics (AISTATS 2016)
% -----------------------------
% Code written by
% Shaobo Han, Duke University
% shaobohan@gmail.com
% 09/29/2015

% Methods
% 1 Ground Truth (Gibbs Sampler)
% 2 Mean-field VB
% 3 Deterministic VGC with Log-normal Margins
% 4 Deterministic VGC with Log-normal Margins and Diagonal Covariance
% 5 Stochastic VGC with Bernstein Polynomials

% We used minFunc package for in deterministic implementation of VGC-LN,
% M. Schmidt. minFunc: unconstrained differentiable multivariate
% optimization in Matlab.
% http://www.cs.ubc.ca/~schmidtm/Software/minFunc.html, 2005.
% To run this code, make sure minFunc is installed

clear all; close all; clc;

 addpath(genpath(pwd))     % Add all sub-directories to the Matlab path

fix. P = 2; % Number of latent variables

trueModel = 'Horseshoe'; fix.y = 0.01;
inferModel = 'MVN'; % Multivariate Gaussian

fix.c = 1;  % unnormalizing constant
%%
opt.k = 10; % Degree/Maximum Degree of Bernstein Polynomials
opt.MaxIter = 50; % Number of SGD stages
opt.NumberZ = 1; % Average Gradients
opt.InnerIter = 500; % Number of iteration
opt.N_mc = 1; % Number of median average in sELBO stepsize search
opt.nsample = 1e6;
figplot.nbins = 60;
BPtype = 'BP';  % Bernstein Polynomials
% BPtype = 'exBP';  % Extended Bernstein Polynomials
% PsiType = 'Normal';   opt.PsiPar(1) = 0;  opt.PsiPar(2) = 1;  % variance
PsiType = 'Exp';    opt.PsiPar = 0.1;
PhiType = 'Normal';  opt.PhiPar(1) = 0; opt.PhiPar(2) =1;  % variance

VGmethod = 'Numeric';

% Learning rate
opt.LearnRate.Mu = 1e-3; opt.LearnRate.C =1e-3;
opt.LearnRate.W = 0.5*1e-3; opt.LearnRate.dec = 0.95; % decreasing base learning rate

switch BPtype
    case 'BP'
        opt.D = opt.k;
    case 'exBP'
        opt.D = opt.k*(opt.k+1)/2; % # of basis functions
end

% Diagonal constraint on Upsilon
opt.diagUpsilon = 0;
opt.adaptivePhi = 0;
opt.normalize = 1;

%% Initialization
ini.Mu = opt.PhiPar(1).*ones(fix.P,1);
% ini.C = sqrt(opt.PhiPar(2)).*eye(fix.P);
ini.C  = randchol(fix.P);
% ini.w = randBPw(fix.P, opt.D, 1, 1);
ini.w = 1./opt.D.*ones(fix.P, opt.D);
% Median Outlier Removal
% Use modified Z score to determine and remove outliers
opt.OutlierTol = 100; % Threshold for online outlier detection
opt.WinSize = 20; % Size of the window
ini.WinSet = ini_WinSet(trueModel, PhiType, PsiType, BPtype, fix, opt, ini);
opt.Wthreshold = 1e12;

%% 1. Gibbs Sampler
Y_mc = Horseshoe_Gibbs(fix.y, opt.nsample);
MCMC = hist1D2D(log(Y_mc), figplot.nbins);
%% 2. MFVB
iterMax = 500; tol = 1e-6;
[Y_vb,ELBO_vb]=Horseshoe_Mfvb(fix.y, opt.nsample, iterMax, tol);
MFVB = hist1D2D(log(Y_vb), figplot.nbins);

%% 3. VC(D): Log-normal
% -------------------full covariance (+C)-------------------
nsample = 30000; nbins =50;
post.m = ini.Mu;
post.c = ini.C;
c1 = 1+0.5*log(2*pi)-2*log(gamma(0.5));
c0 = -0.5*log(2*pi)-2*log(gamma(0.5));

fix.c1 = c1; fix.x = fix.y; fix.c0 = c0;
fix.post.dim = fix.P;

options  =  [];
% options.display  =  'full';
options.display  =  'none';
options.maxFunEvals  =  2000;
options.Method  =  'qnewton';
fix.Flagdiag = 0;  % full covariance
[Y_vcdlnc,ELBO_vcdlnc,postoptimized] = Horseshoe_VCD_LN(post, options,nsample, fix);
VCDLNC = hist1D2D(log(Y_vcdlnc), nbins);
%% 4. VC(D): Log-normal + diagonal
% -------------------diag covariance (+I)-------------------
fix.Flagdiag = 1; % diagonal covariance
[Y_vcdlni,ELBO_vcdlni] = Horseshoe_VCD_LN(post, options,nsample, fix);
VCDLNI = hist1D2D(log(Y_vcdlni), nbins);

%% 5.(S)VGC-BP
fix.scale = 1; % The scale of covariance
inferModel2 = 'MVN'; % Multivariate Normal
% ini.Mu = VG.mu;
% ini.C = VG.C;
opt.updateMuC=1; opt.updateW=1;
[ELBO1, par]   = vgcbp_MuCW(trueModel, inferModel2, PhiType, PsiType, BPtype,VGmethod, fix, ini, opt);
figplot.Y_svc2 = sampleGC(PhiType, PsiType, BPtype, opt, par);
figplot.Y_svc2(sum(~isfinite(log(figplot.Y_svc2)),2)~=0,:) = [];
VGCBPall = hist1D2D(log(figplot.Y_svc2), figplot.nbins);

%% Plot
contourvector = [0.1, 0.25,0.5, 0.75, 0.9];
T1 =MCMC.count';  F1 =T1./sum(sum(T1));  plevs1=contourvector.*max(F1(:));
T2 =MFVB.count';  F2 =T2./sum(sum(T2));  plevs2=contourvector.*max(F2(:));
T3 =VCDLNC.count';  F3 =T3./sum(sum(T3));  plevs3=contourvector.*max(F3(:));
T4 =VCDLNI.count';  F4 =T4./sum(sum(T4));  plevs4=contourvector.*max(F4(:));
T5 =VGCBPall.count';  F5 =T5./sum(sum(T5));  plevs5=contourvector.*max(F5(:));

figure(1)
fs = 16; ncontour = 4; ls = 1.5; ms = 2.5;
main=subplot(4,4,[5,6,7,9,10,11,13,14,15]);
set(gca, 'fontsize',fs)
contour(MCMC.Seqy,MCMC.Seqx,MCMC.count'./sum(sum(MCMC.count)),ncontour, 'Color', 'k', 'LineWidth', 1.2*ls);
hold on
contour(MFVB.Seqy,MFVB.Seqx,MFVB.count'./sum(sum(MFVB.count)),ncontour, '--', 'Color', 'r', 'LineWidth', 2*ls);
hold on
contour(VCDLNC.Seqy,VCDLNC.Seqx,VCDLNC.count'./sum(sum(VCDLNC.count)), ncontour, 'Color', 'm', 'LineWidth', 2*ls);
hold on
contour(VCDLNI.Seqy,VCDLNI.Seqx,VCDLNI.count'./sum(sum(VCDLNI.count)),ncontour, 'Color', 'c', 'LineWidth', 2*ls);
hold on
contour(VGCBPall.Seqy,VGCBPall.Seqx,VGCBPall.count'./sum(sum(VGCBPall.count)),ncontour, 'Color', 'b', 'LineWidth', 1.5*ls);
hold off
grid on
grid minor
axis([-15,5,-15,5])
legend('Gibbs Sampler', 'MFVB', 'VGCLN-full', 'VGCLN-diag', 'VGC-BP-full')

ylabel('log(\gamma)', 'fontsize',fs)
xlabel('log(\tau)', 'fontsize',fs)
pion=subplot(4,4,[8,12,16]); %right plot (rotated)
set(gca, 'fontsize',fs)
plot(MCMC.x_x1,MCMC.f_x1/trapz(MCMC.x_x1,MCMC.f_x1),  'Color', 'k','linewidth', 1.2*ls, 'Markersize', 4*ms);
hold on
plot(MFVB.x_x1,MFVB.f_x1/trapz(MFVB.x_x1,MFVB.f_x1),'--', 'Color', 'r','linewidth', 2.5*ls, 'Markersize', 4*ms);
hold on
plot(VCDLNC.x_x1,VCDLNC.f_x1/trapz(VCDLNC.x_x1,VCDLNC.f_x1), '-o',  'Color', 'm','linewidth', 1.2*ls, 'Markersize', 2*ms);
hold on
plot(VCDLNI.x_x1,VCDLNI.f_x1/trapz(VCDLNI.x_x1,VCDLNI.f_x1), '-s', 'Color', 'c','linewidth', 1*ls, 'Markersize', 1.2*ms);
hold on
plot(VGCBPall.x_x1,VGCBPall.f_x1/trapz(VGCBPall.x_x1,VGCBPall.f_x1), '-x',  'Color', 'b','linewidth', 1*ls, 'Markersize', 4*ms);
hold off
grid on
grid minor
axis([-15,5,...
    0, 1.05.*max([MCMC.f_x1/trapz(MCMC.x_x1,MCMC.f_x1),...
    MFVB.f_x1/trapz(MFVB.x_x1,MFVB.f_x1),...
    VCDLNC.f_x1/trapz(VCDLNC.x_x1,VCDLNC.f_x1),...
    VCDLNI.f_x1/trapz(VCDLNI.x_x1,VCDLNI.f_x1),])])
legend('Gibbs Sampler', 'MFVB', 'VGC-LN-full', 'VGC-LN-diag', 'VGC-BP-full')
view(90, 270)

poz=subplot(4,4,1:3); %upper plot
set(gca, 'fontsize',fs)
plot(MCMC.x_x2,MCMC.f_x2/trapz(MCMC.x_x2,MCMC.f_x2), 'Color', 'k',  'linewidth', 1.2*ls, 'Markersize', 4*ms);
hold on
plot(MFVB.x_x2,MFVB.f_x2/trapz(MFVB.x_x2,MFVB.f_x2), '--','Color', 'r',  'linewidth', 2.5*ls, 'Markersize', 4*ms);
hold on
plot(VCDLNC.x_x2,VCDLNC.f_x2/trapz(VCDLNC.x_x2,VCDLNC.f_x2),'-o', 'Color', 'm',  'linewidth', 1.2*ls, 'Markersize', 2*ms);
hold on
plot(VCDLNI.x_x2,VCDLNI.f_x2/trapz(VCDLNI.x_x2,VCDLNI.f_x2), '-s', 'Color', 'c','linewidth', 1*ls, 'Markersize', 1.2*ms);
hold on
plot(VGCBPall.x_x2,VGCBPall.f_x2/trapz(VGCBPall.x_x2,VGCBPall.f_x2), '-x', 'Color', 'b',  'linewidth', 1*ls, 'Markersize', 4*ms);
hold off
grid on
grid minor
axis([-15,5,...
    0, 1.05.*max([MCMC.f_x2/trapz(MCMC.x_x2,MCMC.f_x2),...
    MFVB.f_x2/trapz(MFVB.x_x2,MFVB.f_x2),...
    VCDLNC.f_x2/trapz(VCDLNC.x_x2,VCDLNC.f_x2),...
    VCDLNI.f_x2/trapz(VCDLNI.x_x2,VCDLNI.f_x2) ])])
pos1=get(poz,'Position'); pos2=get(main,'Position'); pos3=get(pion,'Position');
pos1(3) = pos2(3); %width for the upper plot
set(poz,'Position',pos1)
pos3(4) = pos2(4); %height for the right plot
set(pion,'Position',pos3)
