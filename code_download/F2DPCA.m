function [ P ] = F2DPCA( A_matrix,trainNum,row,column,Num_projection)
%% F2DPCA
%% Wang Q. and Gao Q., et al, Two-dimensional PCA with F-norm minimization, AAAI, 2017:2718-2724. 


%% inputs:
%%  A_matrix -- row*column*trainNum data matrix, row is the data dimension, 
%%  and column is the number of data vectors.
%%  trainNum -- the number of image data.
%%  Num_projection -- the dimension of project vector.

    gama = 1e-5;
    
   W = eye( column,Num_projection );
    Dd = zeros(1,trainNum);
    ZZ = zeros( column,column );
J2 = 0;

    for i = 1 : trainNum


%         Dd(i) = 1/norm(U,'fro') ;           %%weighted2
  %         Dd(i) = 1/sqrt(trace(W.'*A_matrix(:,:,i).'*A_matrix(:,:,i) * W)) ;      %%weighted2
%             Dd(i) = 1/trace(W.'*A_matrix(:,:,i).'*A_matrix(:,:,i) * W);           %%weighted1

             T = A_matrix(:,:,i) - A_matrix(:,:,i) * W*W.';       %% D5
%             Dd(i) = 1/(norm(T,'fro'))*(norm(U,'fro'));             %%  D5
 
%             Dd(i) = 1/(norm(T,'fro')+gama);             %%  D1  (AAAI F-2DPCA)
%              Dd(i) = 1/(norm(T,'fro').^0.5+gama);             %%  Dp  
               Dd(i) = 1/((norm(T,'fro')).^2+gama);             %%  D3(CVPR-workshop)
%                  Dd(i) = 1/((norm(T,'fro')).^3+gama);             %%  D4
         ZZ = ZZ + A_matrix(:,:,i).'*  Dd(i) * A_matrix(:,:,i); clear D;
         J2 = J2 + (norm(T,'fro'));
    end  
     J0 =trace( W.' * ZZ * W );
     JV = J0;
     JJ = J2;
     JD=1;
     
% while  J_D > Precision | numberIteration <= numberIterationMax
    for Ite = 1 :50
%         [H,Sigma,V] = svd( ZZ * W );
%         W = H * eye( column,Num_projection ) * V.'; clear H; clear Sigma; clear V;
%%      H = W.' * ZZ
        [H,Sigma,V] = svd( W.' * ZZ );  
        W = V * eye( column,Num_projection ) * H.'; clear H; clear Sigma; clear V;
        
     J2 = 0;
    Dd = zeros(1,trainNum);
    ZZ = zeros( column,column );
    for i = 1 : trainNum
        U = A_matrix(:,:,i) * W;
        
%         Dd(i) = 1/norm(U,'fro') ;           %%weighted2
  %         Dd(i) = 1/sqrt(trace(W.'*A_matrix(:,:,i).'*A_matrix(:,:,i) * W)) ;      %%weighted2
%             Dd(i) = 1/trace(W.'*A_matrix(:,:,i).'*A_matrix(:,:,i) * W);           %%weighted1
         
            T = A_matrix(:,:,i) - A_matrix(:,:,i) * W*W.';       %% D5
%              Dd(i) = 1/(norm(T,'fro'))*(norm(U,'fro'));             %%  D5

              Dd(i) = 1/(norm(T,'fro')+gama);             %%  D1 (AAAI F-2DPCA)
%              Dd(i) = 1/(norm(T,'fro').^0.5+gama);             %%  Dp  
%               Dd(i) = 1/((norm(T,'fro')).^2+gama);             %%  D3 (CVPR-workshop)
%                  Dd(i) = 1/((norm(T,'fro')).^3+gama);             %%  D4
         ZZ = ZZ + A_matrix(:,:,i).'*  Dd(i) * A_matrix(:,:,i); clear D;
                J2 = J2 + (norm(T,'fro'));
    end  
    
      J1 =trace( W.' * ZZ * W );
        JV=[JV J1];
        JD = J1-J0;
        J0=J1;  clear J1;

  JJ=[JJ J2];
  
    end
    figure(1)
    plot(JV);
P = W;