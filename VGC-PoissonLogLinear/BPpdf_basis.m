% Bernstein polynomial PDF basis function
% Beta densities
% u: evaluating point, could be a p by 1 vector 
% k: degree of BPs
% opt: BP/exBP, standard or extended BP

% V: a p by d matrix storing basis funcitons 

% Shaobo Han
% 08/08/2015

function V = BPpdf_basis(u, k, opt)
if nargin<3,
    opt = 'BP'; % Standard BP by default
end
p = length(u); % # of variables
switch opt
    case 'BP' % Standard BP
        d = k; % # of basis functions
        a_seq = 1:k; b_seq = k-a_seq+1;
    case 'exBP' % extended BP
        d = k*(k+1)/2; % # of basis functions
        a_seq = zeros(1, d); b_seq = zeros(1, d);
        tmp_idx = cumsum(1:k);
        for j = 1:k
            tmp2=tmp_idx(j); tmp1=tmp_idx(j)-j+1;
            tmp_a = 1:j; tmp_b = j-tmp_a+1;
            a_seq(1, tmp1:tmp2) = tmp_a; b_seq(1, tmp1:tmp2) = tmp_b;
        end
end
A_seq = repmat(a_seq, p,1); B_seq = repmat(b_seq, p,1);
U = repmat(u, 1, d); % p by d matrix
V = betapdf(U , A_seq, B_seq); % p by d matrix
end
