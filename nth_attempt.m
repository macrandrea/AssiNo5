clear
clc
maxit=30;
tol=1e-6;
r=ones(maxit,8);
x=ones(maxit,2);
l=ones(maxit,3);
z=ones(maxit,3);
mu=1;
x(1,1:2)=[1.9,1.44];
l(1,1:3)=[2,3,4];
z(1,1:3)=[1,3,4];
Z=[z(1,1),0,0;0,z(1,2),0;0,0,z(1,3)];
L=[l(1,1),0,0;0,l(1,2),0;0,0,l(1,3)];
I=eye(3);
A=[-1,-1;1,0;0,1];
H=[2 , 0;0 , 2,];
passoX=1;
passoL=1;
passoZ=1;
for i=1:maxit
    r(i,1)=2*(x(i,1)-4)-l(i,1)+l(i,2);%deriv Lagr di x1    
    r(i,2)=2*x(i,2)-l(i,1)+l(i,3);    %deriv Lagr di x2     
    r(i,3)=2-x(i,1)-x(i,2)-z(i,1);    %deriv Lagr di l1  
    r(i,4)=x(i,1)-z(i,2);         %deriv Lagr di l2 
    r(i,5)=x(i,2)-z(i,3);         %deriv Lagr di l3 
    r(i,6)=z(i,1)*l(i,1)-mu;      %deriv Lagr di z1  
    r(i,7)=z(i,2)*l(i,2)-mu;      %deriv Lagr di z2     
    r(i,8)=z(i,3)*l(i,3)-mu;      %deriv Lagr di z3
    deltaX=(H+A'*(Z\L)*A)\(-r(i,1:2)'-A'*(inv(Z)*(r(i,3:5)'+inv(L)*r(i,6:8)')));
    deltaL=(L\Z)\(-A*deltaX-r(i,3:5)'+L\r(i,6:8)');
    deltaZ=-L\(r(i,6:8)'+Z*deltaL);
    x(i+1,:)=x(i,:)+0.01*deltaX';
    l(i+1,:)=l(i,:)+0.01*deltaL';
    z(i+1,:)=z(i,:)+0.01*deltaZ';
    if(norm(r(i,:))<tol)
        disp('break')
        break;
    else
        mu=mu/2;
    end
end