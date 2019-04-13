% To implement the example in Sec. 3.5.D for classifying handwritten digits from MNIST database.
% Input:
% X: Ten classes of input data, each is of size 784 by nt (here we use nt = 500).
% n: the actual size of each data matrix is 784 by n. Hence n is limited to n <= nt.
% K: the number of "eigen-digits" to be used in PCA. Typically K << m = 784.
% Te: testing data (here we use Te28 of size 784 by 10000).
% Lte: Labels for the testing data.
% Output:
% Js: the numerical values of testing digits obtained by PCA-based classification.
% er: classification error.
% Written by W.-S. Lu, University of Victoria. Last modified: March 9, 2015.
% Example: load X500; load Te28; load Lte28;
% [Js,er] = class_pca_digits(X500,500,24,Te28,Lte28);
function [Js,er] = class_pca_digits(X,n,K,Te,Lte)
U = [];
Xb = [];
t = size(X,2);
nt = t/10;
for i = 1:10,
    Tw = X(:,((i-1)*nt+1):(i*nt));
    Ti = Tw(:,1:n);
    xbi = mean(Ti')';
    Xh = Ti - xbi*ones(1,n);
    [u,s,v] = svd(Xh*Xh');
    U = [U u(:,1:K)];
    Xb = [Xb xbi];
end
Lte = Lte(:);
M = length(Lte);
Js = zeros(M,1);
for m = 1:M,
    xw = Te(:,m);
    e = zeros(1,10);
    for j = 1:10;
        Uj = U(:,((j-1)*K+1):(j*K));
        xbj = Xb(:,j);
        zj = Uj'*(xw-xbj);
        xjh = Uj*zj + xbj;
        e(j) = norm(xw-xjh);
    end
    [emin,ind] = min(e);
    Js(m) = ind - 1;
end
E = (Lte ~= Js);
er = sum(E)/M