function [ELBO, par]  = vgcbp_MuCW(trueModel, inferModel, PhiType, PsiType, BPtype, VGmethod, fix, ini, opt)
% VGmethod = 'Analytic' | 'Numeric';

% Parameters Unpack
InnerIter = opt.InnerIter;  MaxIter = opt. MaxIter;
r1 = opt.LearnRate.Mu;  r2 = opt.LearnRate.C;
r3 = opt.LearnRate.W; dec = opt.LearnRate.dec;

if opt.diagUpsilon ==1
    display(['[VGC(S): diag], BPtype = ', BPtype,' PsiType = ', PsiType,' PhiType = ', PhiType])
else
    display(['[VGC(S): full], BPtype = ', BPtype,' PsiType = ', PsiType,' PhiType = ', PhiType])
end
ELBO = zeros(InnerIter, MaxIter);
RMSE = zeros(InnerIter, MaxIter);
for iter = 1: MaxIter
    tic;
    for inner = 1:InnerIter
        % First Iteration
        if (iter == 1)&&(inner==1)
            par.Mu_old = ini.Mu; par.C_old  =  ini.C; par.W_old =  ini.w;
            WinSet = ini.WinSet;
            par.Mu = par.Mu_old; par.C = par.C_old; par.W = par.W_old;
        else
            par.Mu_old = par.Mu; par.C_old = par.C; par.W_old = par.W;
        end
        
        if opt.updateMuC==1
        %%  Draw Samples and Calculate (Mu, C)
        [Delta_Mu, Delta_C, WinSet] = Grad_MuC(VGmethod, trueModel, inferModel, PhiType, PsiType, BPtype, fix, opt, par,WinSet);
               if opt.normalize ==1
        if norm(Delta_Mu,'fro')<=1e-5
            Delta_Mu = zeros(P,1);
        else
            Delta_Mu = Delta_Mu./norm(Delta_Mu,'fro');
        end
        if norm(Delta_C,'fro')<=1e-5
            Delta_C = zeros(P,P);
        else
            Delta_C = Delta_C./norm(Delta_C,'fro');
        end
               end
        %% Update Mu
        par.Mu = par.Mu_old + r1.*Delta_Mu;
        if opt.adaptivePhi == 1
            par.Mu = zeros(size( par.Mu));
        end
        %% Update Mu
        par.C = par.C_old + r2.*Delta_C;
        else
        par.Mu = par.Mu_old;  par.C = par.C_old;   
        end
        if opt.updateW==1
        %%  Draw Samples and Calculate Delta_W
        [Delta_w, WinSet] = Grad_W(trueModel, PhiType, PsiType, BPtype, fix, opt, par,WinSet);
         if opt.normalize ==1
        % Normalize
        if norm(Delta_w,'fro')<=1e-5
            Delta_w = zeros(size(Delta_w));
        else
            Delta_w = diag(1./sqrt(sum(abs(Delta_w).^2,2)))*Delta_w;
        end
         end
        %% Update W
        tmpNorm = sqrt(sum(abs(Delta_w).^2,2));
        UnstableIndex = (tmpNorm>opt.Wthreshold); 
        tmpsum = sum(UnstableIndex);
        if tmpsum ~=0
            Delta_w(UnstableIndex,:) = zeros(tmpsum, opt.D);
            display('Warning: Unstable Delta W Omitted!')
            tmpNorm(UnstableIndex)
        end       
        W0 = par.W_old + r3.*Delta_w;
        % Project Gradient Descent
        par.W = SimplxProj(W0);
        if sum(abs(sum(par.W,2)-1)>1e-5*ones(size(sum(par.W,2))))>0
            display('Error: Weight not Sum up to 1!')
            sum(par.W,2)
            par.W = par.W_old;
        end     
        else
           par.W = par.W_old;    
        end
        ELBO(inner, iter) = Cal_mcELBO(trueModel, inferModel, PhiType, PsiType, BPtype,  fix, opt, par);
        RMSE(inner, iter) = norm(diag(corrcov(par.C *par.C'),-1)- diag(fix.trueUpsilon,-1), 'fro')./norm(diag(fix.trueUpsilon,-1), 'fro');
    end
    r1 = r1*dec;     r2 = r2*dec;     r3 = r3*dec;
    t1= toc;
    display(['ELBO = ', num2str(median(ELBO(:,iter))),  ', Iter: ' num2str(iter), '/' num2str(MaxIter),' Elapsed: ', num2str(t1) ' sec']);
end
par.WinSet = WinSet;
par.RMSE = RMSE; 
end
