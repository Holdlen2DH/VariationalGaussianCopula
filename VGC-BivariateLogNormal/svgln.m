function output = svgln(truemodel, infermodel, stageMax, iterMax, method, fix)
% method = 'Analytic' | 'Numeric';

P = fix.P; rho = fix.rho;
dec = fix.dec; scale = fix.scale;
% Initialize
tmpMu = zeros(P, 1);
tmpC = sqrt(scale).*0.1.*eye(P);

ELBO = zeros(stageMax, iterMax);
RMSE1 = zeros(stageMax, iterMax);
RMSE2 = zeros(stageMax, iterMax);

for stage = 1: stageMax
    for iter = 1:iterMax
        if (iter ==1) && (stage==1)
            Mu = tmpMu; C = tmpC;
        end
        tmpw = randn(P,1);
        tmpz = Mu + C*tmpw;
        switch method
            case 'Analytic'
                tmpx = exp(tmpz);
                g = derivatives(tmpx, truemodel, fix).*exp(tmpz)+ones(2,1);
                dmu = g;
                dC = tril(g*tmpw') + diag(1./diag(C));
            case 'Numeric'
                tmpx = exp(tmpz);
                g1 = derivatives(tmpx, truemodel, fix).*exp(tmpz)+ones(2,1);
                par.Mu = Mu; par.Sigma = C*C';
                g2 = logqpost_derivatives(tmpz, infermodel, par);
                dmu = g1-g2;
                dC = tril(dmu*tmpw');
        end
        Mu = Mu +  rho.*dmu;
        C = C +  rho.*dC;
        par.Mu = Mu; par.Sigma = C*C'; par.C = C;
        tmpw = randn(P,1);
        tmpz = Mu + C*tmpw;
        tmpx = exp(tmpz);
        sELBO = logmodel(tmpx, truemodel, fix)+sum(tmpz);
        logq = logqpost(tmpz, infermodel, par);
        ELBO(stage, iter) = sELBO-logq;
        RMSE1(stage, iter) = norm(Mu- fix.trueMu, 'fro')./norm(fix.trueMu, 'fro');
        RMSE2(stage, iter) = norm(diag(corrcov(C*C'),-1)-  diag(fix.trueUpsilon,-1), 'fro')./norm(diag(fix.trueUpsilon,-1), 'fro');
    end
    rho = rho*dec;
end
fELBO = median(ELBO(end, :));
output.RMSEmu = RMSE1;
output.RMSECTC = RMSE2;
output.fsELBO = fELBO;
output.sELBO = ELBO;
output.mu = Mu;
output.C = C;
end