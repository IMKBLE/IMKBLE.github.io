function [ P,M ] = OM_F2DPCA( A_matrix,trainNum,row,column,Num_projection)
%% This routine solves the following problem, which has been also presented 
%% in the paper entitled "Optimal mean two-dimensional principal  
%% component analysis with F-norm minimization".
%% Wang Q., Gao Q., et al. Optimal mean two-dimensional principal component analysis with F-norm minimization, Pattern Recognition, 2017, 68(C): 286-294. 


%% inputs:
%%  A_matrix -- row*column*trainNum data matrix, row is the data dimension, 
%%  and column is the number of data vectors.
%%  trainNum -- the number of image data.
%%  Num_projection -- the dimension of project vector.


W = eye( column,Num_projection );
D = ones( 1,trainNum );
gama = 1e-5;
numberIteration = 0;
numberIterationMax = 80;
Precision = 1e-8;
sumD = sum(D);
sumDA = zeros( row,column );
sumT = zeros( column,column );

for i = 1 : trainNum
    U = D(i)* A_matrix(:,:,i) ;
    sumDA = sumDA +  U;
end
M = sumDA/sumD; clear sumDA U sumD;

for i = 1 : trainNum
    U = A_matrix(:,:,i)-M;
    T = U.'* D(i) *U;
    sumT = sumT  +  T;
end
clear U T;

J0 = trace(W.'* sumT * W); 
JV1 = J0;
J_value = 1;
 for Ite = 1 :30
%     numberIteration = numberIteration+1;
%% update W
    [H,Sigma,V] = svd(sumT * W );
    W = H * eye( column,Num_projection ) * V.'; clear H; clear Sigma; clear V;
%% update D    
    for i = 1 : trainNum
        U = A_matrix(:,:,i)-M ;
        E = U - U * W * W';
        NE = norm(E,'fro');
        D(i) = 1/(NE+gama);
    end
    
    sumD = sum(D);
    sumDA = zeros( row,column );
    sumT = zeros( column,column );
%% update M        
    for i = 1 : trainNum
        U = D(i)* A_matrix(:,:,i) ;
        sumDA = sumDA + U;
    end
    M = sumDA/sumD;clear sumDA U sumD;
    
    for i = 1 : trainNum
        U = A_matrix(:,:,i)-M ;
        T = U.'* D(i) *U;
        sumT = sumT +  T;
    end
    
    J1 = trace(W.'* sumT * W); clear U T;
    JV1=[JV1 J1];
    J_value = J1-J0;
    J0=J1;
    clear J1;   
end
figure(1)
plot(JV1);
P = W;