function [W,JJ ]= L2_WeightedPCAMinSUM( A,dim)
%% weightedPCA
%% Wang Q, Gao Q, Gao X, et al. L2,p-norm based PCA for Image Recognition.
%% IEEE Transactions on Image Processing, 2017, PP(99):1-1.
% A :data r*c
% dim£ºdimension number

[r,c]= size(A)
 vvv = rand(r,dim);
% initial
% load vvv
W=orth(vvv);

numberIteration = 1;
numberIterationMax = 50;
I = eye(r,dim);

J_value = 1;
J0 = 0;
JJ = [];
lambda = 0.0001;
Precision = 1e-3;

%iteration
% while  J_value > Precision && numberIteration < numberIterationMax
while  numberIteration <= numberIterationMax
    d = zeros(1,c);
    d3 = zeros(1,c);
  for i = 1:c
      d1 = norm((A(:,i)-W*W'*A(:,i)),2);
      d2 = norm((W'*A(:,i)),2);
%     tr = trace(W'*A(:,i)*A(:,i)'*W);
%       d(i) = 1/ (d1^2 + lambda);   %% p = 0
%       d(i) = 1/ (d1^1.5 + lambda); %% p = 0.5
      d(i) = 1/ (d1 + lambda); %% p = 1

  end

  D = diag(d);
  H = A*D*A'*W;
  
J1 = trace(W'*H);
JJ = [JJ,J1];

J_value = J1 - J0;
J0 = J1;

[U,S,V] = svd(H);

W = U*I*V';
numberIteration = numberIteration +1;
end

figure(1)
plot(JJ)


