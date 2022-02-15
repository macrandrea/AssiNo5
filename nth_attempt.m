clear
clc
maxit=30;
tol=1e-6;
r=ones(maxit,8);
x=ones(maxit,2);
l=ones(maxit,3);
z=ones(maxit,3);
mu=10;
x(1,1:2)=[0.89,0.01];
l(1,1:3)=[4,0.01,4];
z(1,1:3)=[mu/4,0.01,mu/4];
passoX=1;
passoL=1;
passoZ=1;
for i=1:maxit
    while(norm(r(i,:))>tol)
    r(i,1)=2*(x(i,1)-4)-l(i,1)+l(i,2);%deriv Lagr di x1    
    r(i,2)=2*x(i,2)-l(i,1)+l(i,3);    %deriv Lagr di x2     
    r(i,3)=2-x(i,1)-x(i,2)-z(i,1);    %deriv Lagr di l1  
    r(i,4)=x(i,1)-z(i,2);         %deriv Lagr di l2 
    r(i,5)=x(i,2)-z(i,3);         %deriv Lagr di l3 
    r(i,6)=z(i,1)*l(i,1)-mu;      %deriv Lagr di z1  
    r(i,7)=z(i,2)*l(i,2)-mu;      %deriv Lagr di z2     
    r(i,8)=z(i,3)*l(i,3)-mu;      %deriv Lagr di z3
    menoR=-r(i,:);
    J=[2 , 0 , -1, 1   , 0 , 0 ,   0 , 0;
       0 , 2 , -1, 0   , 1 , 0 ,   0 , 0;
       1 , 1, 0  , 0   , 0 , -1,   0 , 0;
       1 , 0 , 0 , 0   , 0 , 0 ,   -1, 0;
       0 , 1 , 0 , 0   , 0 , 0 ,   0 ,-1;
       0 , 0 , z(i,1), 0 , 0 , l(i,1), 0 , 0;
       0 , 0 , 0 , z(i,2), 0 , 0 , l(i,2), 0;
       0 , 0 , 0 , 0 , z(i,3), 0 , 0 ,l(i,3);
       ];
    dd=1;
    while (dd>tol)
    delta=J\menoR';%congj gradient with tol on relative res
    res=J*r(i,:)'-delta;
    dd=norm(res)-0.9*l(i,:)*z(i,:)';
    end
    x(i+1,:)=x(i,:)+1*delta(1:2)';
    l(i+1,:)=l(i,:)+1*delta(3:5)';
    z(i+1,:)=z(i,:)+1*delta(6:8)';
    %if(norm(r(i,:))<tol)
    %    disp('break')
    %    break;
    %else
        mu=mu/2;
    end
end