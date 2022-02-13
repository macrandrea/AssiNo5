clear, clc
%initialization
maxit=100;
tol=1e-6;
mu=20;
r=ones(8,maxit)';
mat=ones(8,maxit)';
x1=mat(1,1);
x2=mat(1,2);
l1=mat(1,3);
l2=mat(1,4);
l3=mat(1,5);
z1=mat(1,6);
z2=mat(1,7);
z3=mat(1,8);
delta=zeros(8,maxit)';
f=zeros(maxit,1);
f(1)=(mat(1,1)-4)^2-mat(1,2)^2;
alpha=ones(maxit,8);
    r(1,1)=2*(x1-4)-l1+l2;%deriv Lagr di x1    
    r(1,2)=2*x2-l1+l3;    %deriv Lagr di x2     
    r(1,3)=2-x1-x2-z1;    %deriv Lagr di l1  
    r(1,4)=x1-z2;         %deriv Lagr di l2 
    r(1,5)=x2-z3;         %deriv Lagr di l3 
    r(1,6)=z1*l1-mu;      %deriv Lagr di z1  
    r(1,7)=z2*l2-mu;      %deriv Lagr di z2     
    r(1,8)=z3*l3-mu;      %deriv Lagr di z3  

for i=1:maxit
    
    J=[2 , 0 , -1, 1   , 0 , 0 ,           0 , 0;
       0 , 2 , -1, 0   , 1 , 0 ,           0 , 0;
       1 , 1, 0  , 0   , 0 , -1,           0 , 0;
       1 , 0 , 0 , 0   , 0 , 0 ,           -1, 0;
       0 , 1 , 0 , 0   , 0 , 0 ,           0 ,-1;
       0 , 0 , mat(i,6), 0 , 0 , mat(i,3), 0 , 0;
       0 , 0 , 0 , mat(i,7), 0 , 0 , mat(i,4), 0;
       0 , 0 , 0 , 0 , mat(i,8), 0 , 0 ,mat(i,5);
       ];
    delta(i,:)=cgs(J,r(i,:)',tol); %conjugate gradient with tolerance tol
    mat(i+1,:)=mat(i,:)-alpha(i,:)*delta(i,:)';
    
    
        r(i,1)=2*(mat(i+1,1)-4)-mat(i+1,3)+mat(i+1,4);%deriv Lagr di x1    
        r(i,2)=2*mat(i+1,2)-mat(i+1,3)+mat(i+1,5);    %deriv Lagr di x2     
        r(i,3)=2-mat(i+1,1)-mat(i+1,2)-mat(i+1,6);    %deriv Lagr di l1  
        r(i,4)=mat(i+1,1)-mat(i+1,7);         %deriv Lagr di l2 
        r(i,5)=mat(i+1,2)-mat(i+1,8);         %deriv Lagr di l3 
        r(i,6)=mat(i+1,6)*mat(i+1,3)-mu;      %deriv Lagr di z1  
        r(i,7)=mat(i+1,7)*mat(i+1,4)-mu;      %deriv Lagr di z2     
        r(i,8)=mat(i+1,8)*mat(i+1,5)-mu;      %deriv Lagr di z3  
    if (norm(r)>1e-6)
        break;
    end
end