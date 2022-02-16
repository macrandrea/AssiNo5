clear
clc
maxit=50;
tol=1e-6;
r=ones(maxit,8);
x=ones(maxit,2);
l=ones(maxit,3);
z=ones(maxit+1,3);
mu=20;
x(1,1:2)=[1.9,0.44];
l(1,1:3)=[2,3,4];
z(1,1:3)=[1,3,4];
I=eye(3);
A=[-1,-1;1,0;0,1];
H=[2 , 0;0 , 2,];
passoX=[0.001,0.01]';
e=[1,1,1];
f(1)=(x(1,1)-4)^2-x(1,2)^2+mu*e*log(z(1,1:3))';
for i=1:maxit
    r(i,1)=2*(x(i,1)-4)-l(i,1)+l(i,2);%deriv Lagr di x1    
    r(i,2)=2*x(i,2)-l(i,1)+l(i,3);    %deriv Lagr di x2     
    r(i,3)=2-x(i,1)-x(i,2)-z(i,1);    %deriv Lagr di l1  
    r(i,4)=x(i,1)-z(i,2);         %deriv Lagr di l2 
    r(i,5)=x(i,2)-z(i,3);         %deriv Lagr di l3 
    r(i,6)=z(i,1)*l(i,1)-mu;      %deriv Lagr di z1  
    r(i,7)=z(i,2)*l(i,2)-mu;      %deriv Lagr di z2     
    r(i,8)=z(i,3)*l(i,3)-mu;      %deriv Lagr di z3
    Z=[z(i,1),0,0;0,z(i,2),0;0,0,z(i,3)];
    L=[l(i,1),0,0;0,l(i,2),0;0,0,l(i,3)];
    deltaX=(H+A'*(Z\L)*A)\(-r(i,1:2)'-A'*(Z\r(i,3:5)'+L\r(i,6:8)'));
    deltaL=(L\Z)\(-A*deltaX-r(i,3:5)'+L\r(i,6:8)');
    deltaZ=-L\(r(i,6:8)'+Z*deltaL);
    slX=deltaX'*r(i,1:2)';
    %if(norm(r(i,:))>tol)
        x(i+1,:)=x(i,:)+0.01*deltaX';
        l(i+1,:)=l(i,:)+0.01*deltaL';
        z(i+1,:)=z(i,:)+0.01*deltaZ';  
        f(i+1)=(x(i+1,1)-4)^2-(x(i+1,2))^2-mu*e*log(z(i+1,:))';
        if(f(i+1)>f(i))
            mu=mu/2;
        else
            break
        end
    if(norm(r(i,:))>tol)
        mu=mu/2;
    else
       break;
    end
end