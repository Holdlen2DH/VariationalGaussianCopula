function [ELBO, par] = vgcbp_w(trueModel, inferModel,  PhiType, PsiType, BPtype, fix, ini, opt)
% Parameters Unpack
InnerIter = opt.InnerIter;  MaxIter = opt. MaxIter;
r3 = opt.LearnRate.W;  dec = opt.LearnRate.dec;
if opt.diagUpsilon ==1
    display(['[VGC(S): diag], BPtype = ', BPtype,' PsiType = ', PsiType,' PhiType = ', PhiType])
else
    display(['[VGC(S): full], BPtype = ', BPtype,' PsiType = ', PsiType,' PhiType = ', PhiType])
end
ELBO = zeros(InnerIter, MaxIter);
for iter = 1: MaxIter
    tic;
    for inner = 1:InnerIter
        % First Iteration
        if (iter == 1)&&(inner==1)
            par.Mu_old = ini.Mu; par.C_old  =  ini.C; par.W_old =  ini.w;
            WinSet = ini.WinSet;
        else
            par.Mu_old = par.Mu; par.C_old = par.C; par.W_old = par.W;
        end
        % Do Not Update (Mu, C) Here
        par.Mu = par.Mu_old;  par.C = par.C_old;
        %%  Draw Samples and Calculates
        [Delta_w, WinSet] = Grad_W(trueModel, PhiType, PsiType, BPtype, fix, opt, par,WinSet);
        if opt.normalize ==1
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
            display('Error: weight not sum up to 1!')
            sum(par.W,2)
            par.W = par.W_old;
        end
        ELBO(inner, iter) = Cal_mcELBO(trueModel, inferModel, PhiType, PsiType, BPtype,  fix, opt, par);
    end
    r3 = r3*dec;
    t1= toc;
    display(['ELBO = ', num2str(median(ELBO(:,iter))),  ', Iter: ' num2str(iter), '/' num2str(MaxIter),' Elapsed: ', num2str(t1) ' sec']);
end
end