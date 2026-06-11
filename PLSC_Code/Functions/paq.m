function [P,a,Q]=paq(X,k)
%USAGE [P,a,Q]=paq(X,k);
% k-Singular value decomposition of the I by J matrix X
%if k is not present then k=min{I,J}
%if k is larger larger than min{I,J} then k= min
%if k is larger than K=the actual number of singular values, then k=K
% the eigenvectors and singular-values are ordered
% in decreasing order
% P are the eigenvectors of X'X
% Q are the eigenvector of XX'
% a is the vector of the SINGULAR values
% NOTE that a = sqrt(lambda)
% where lambda are the eigenvalues of both X'X and XX'
% 
[I,J]=size(X);
m=min(I,J); 
if nargin==1, k=m;
        else if k > m, k=m;end;
     end;
flip=0; 
if I < J, X=X';flip=1;end;
[Q,a]=eigen(X'*X);
l= max(size(a)); if k > l, k=l;end;
Q=Q(:,1:k);
a=a(1:k);
a=sqrt(a);
P=X*Q*inv(diag(a));
if flip==1,X=X';
 bidon=Q;Q=P;P=bidon;end