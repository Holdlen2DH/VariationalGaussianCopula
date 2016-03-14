%  psi function: small psi, derivative of Psi
% Shaobo Han
% 08/10/2015
function output = spsi_f4(x, opt)
output = zeros(4,1); 
output(1)= spsi_f(x(1), opt.PoM.PsiPar1, opt.PoM.PsiType1); 
output(2)= spsi_f(x(2), opt.PoM.PsiPar1, opt.PoM.PsiType1); 
output(3)= spsi_f(x(3), opt.PoM.PsiPar1, opt.PoM.PsiType1); 
output(4)= spsi_f(x(4), opt.PoM.PsiPar2, opt.PoM.PsiType2); 
end