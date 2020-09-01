function [W]= Angle_PCA( A,dim)
%% angle PCA
%% Wang Q, Gao Q, Gao X, et al.Angle principal component analysis, IJCAI, 2017:2936-2942.
 [r,c] = size(A);
 vvv = rand(r,dim);
% 初始化
% load vvv
W=orth(vvv);
% d = zeros(1,c);
numberIteration = 1;
numberIterationMax = 20;
% I = eye(r,dim);

J_value = 1;
J0 = 0;
JJ = [];
gama = 1;
Precision = 1e-8;
%迭代更新W
% while  J_value > Precision && numberIteration < numberIterationMax
while  numberIteration < numberIterationMax
%       for numberIteration = 1:numberIterationMax
   
    AQ = sqrt(sum((W'*A).*(W'*A),1));
    BQ = sqrt(sum((A-W*W'*A).*(A-W*W'*A),1));
%     BBQ = BQ;

    id = find (BQ<gama);
    BQ(id) = BQ(id) + gama;

%         BQ = BQ+gama;
    Lamd = AQ./BQ;
    H = A*(diag(1./AQ+Lamd./BQ))*A';
    J1 =sum(Lamd);
    
%   for i = 1:c
%       d1 = norm((A(:,i)-W*W'*A(:,i)),2);
%       d2 = norm((W'*A(:,i)),2);
%       d(i) = d2/(d1+gama);
%       H = H + (1/d2 + d(i)/(d1+gama))*A(:,i)*A(:,i)';      
%       J1 = J1+d2/(d1+gama);
%   end
  
  [eigvector,eigvalue]=eig((H+H')/2);
eigvalue = diag(eigvalue);
% eigIdx = find(eigvalue/max(eigvalue) < 1e-12);
% eigIdx = find(eigvalue < 1e-6);
% eigvalue(eigIdx) = [];
% eigvector(:,eigIdx) = [];
[~, index] = sort((eigvalue),'descend');
% eigvalue = eigvalue(index(1:dim));
W = eigvector(:,index(1:dim));
% W = eigvector;
JJ = [JJ,J1];

J_value = abs(J1 - J0)/J0;
J0 = J1;

% J2 = trace(A*D*A') - J1;
% JJ2 = [JJ2,J2];
numberIteration = numberIteration +1;
end

numberIteration = numberIteration+0
figure(1)
plot(JJ)


% figure(2)
% hold on
% plot(JJ2)

