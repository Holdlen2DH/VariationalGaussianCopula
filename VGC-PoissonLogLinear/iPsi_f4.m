% Invserse Psi function:
% Shaobo Han
% 08/10/2015
function output = iPsi_f4(p, opt)
output = zeros(4,1); 
output(1)= iPsi_f(p(1), opt.PoM.PsiPar1, opt.PoM.PsiType1); 
output(2)= iPsi_f(p(2), opt.PoM.PsiPar1, opt.PoM.PsiType1); 
output(3)= iPsi_f(p(3), opt.PoM.PsiPar1, opt.PoM.PsiType1); 
output(4)= iPsi_f(p(4), opt.PoM.PsiPar2, opt.PoM.PsiType2); 
end