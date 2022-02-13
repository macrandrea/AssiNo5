clear, clc
%initialization
maxit=100;
tol=1e-6;
mu=0.02;
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
alpha=ones(2,1);
for i=1:maxit
    r(1,1)=2*(x1-4)-l1+l2;%deriv Lagr di x1    
    r(1,2)=2*x2-l1+l3;    %deriv Lagr di x2     
    r(1,3)=2-x1-x2-z1;    %deriv Lagr di l1  
    r(1,4)=x1-z2;         %deriv Lagr di l2 
    r(1,5)=x2-z3;         %deriv Lagr di l3 
    r(1,6)=z1*l1-mu;      %deriv Lagr di z1  
    r(1,7)=z2*l2-mu;      %deriv Lagr di z2     
    r(1,8)=z3*l3-mu;      %deriv Lagr di z3  
    J=[2 , 0 , -1, 1   , 0 , 0 ,           0 , 0;
       0 , 2 , -1, 0   , 1 , 0 ,           0 , 0;
       1 , 1, 0  , 0   , 0 , -1,           0 , 0;
       1 , 0 , 0 , 0   , 0 , 0 ,           -1, 0;
       0 , 1 , 0 , 0   , 0 , 0 ,           0 ,-1;
       0 , 0 , mat(i,6), 0 , 0 , mat(i,3), 0 , 0;
       0 , 0 , 0 , mat(i,7), 0 , 0 , mat(i,4), 0;
       0 , 0 , 0 , 0 , mat(i,8), 0 , 0 ,mat(i,5);
       ];
    delta(i,:)=J\r(i,:)';
    %discesa x_1,x_2
    slope=delta(i,1:2)*r(i,1:2)';
    if (slope>0)
        while true
            mat(i+1,1:2)=mat(i,1:2)+alpha'*delta(i,1:2)';
            f(i+1)=(mat(i+1,1)-4)^2-mat(i,2)^2;
            if f(i+1)<f(i)-0.1*alpha*slope
                break
            else
                alpha=alpha/2;
            end
        end
    %aggiungere condizioni di discesa singole e passi singoli
    end
    if(norm(r)<1e-6)
        break;
    end
end