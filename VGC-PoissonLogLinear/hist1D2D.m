function Output=hist1D2D(Matrix, nbins) 
% 1D and 2D histogram
% Matrix N by 2 matrix 
[f_x1,x_x1] = hist(Matrix(:,1), nbins);
[f_x2,x_x2] = hist(Matrix(:,2), nbins);
[count,Dy] = hist3(Matrix,[nbins+2,nbins]);
Seqx= Dy{1,1}; Seqy = Dy{1,2};
Output.f_x1 = f_x1; Output.x_x1 = x_x1; 
Output.f_x2 = f_x2; Output.x_x2 = x_x2; 
Output.count =count'; Output.Seqx =Seqx; Output.Seqy =Seqy; 
end