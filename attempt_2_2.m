maxit=1000;
tol=1e-12;
r=ones(maxit,8);
x=ones(maxit,2);
l=ones(maxit,3);
z=ones(maxit+1,3);
a=ones(maxit,1);
f=ones(maxit+1,1);
mu=1;
% x(1,1:2)=rand(1,2);
% k=rand(1,1)*pi/2;
% x(1,1)=x(1,1)*cos(k);
% x(1,2)=x(1,1)*sin(k);
x(1,:)=[10^-8,10^-8];
l(1,1:3)=[1,1,1];
z(1,1:3)=[1-x(1,1)^2-x(1,2)^2,x(1,1),x(1,2)];
e=[1,1,1];
f(1)=2*(x(1,2))-x(1,1)^2-mu*e*log(z(1,1:3))'-l(1,:)*[1-x(1,1)^2-x(1,2)^2-z(1,1),x(1,1)-z(1,2),x(1,2)-z(1,3)]';
for i=1:maxit
    r(i,1)=-l(i,2)-2*(-l(i,1)+1)*x(i,1);%deriv Lagr di x1    
    r(i,2)=2+2*l(i,1)*x(i,2)-l(i,3);    %deriv Lagr di x2     
    r(i,3)=-1+x(i,1)^2+x(i,2)^2+z(i,1);    %deriv Lagr di l1  
    r(i,4)=-x(i,1)+z(i,2);         %deriv Lagr di l2 
    r(i,5)=-x(i,2)+z(i,3);         %deriv Lagr di l3 
    r(i,6)=z(i,1)*l(i,1)-mu;      %deriv Lagr di z1  
    r(i,7)=z(i,2)*l(i,2)-mu;      %deriv Lagr di z2     
    r(i,8)=z(i,3)*l(i,3)-mu;      %deriv Lagr di z3
    J=[-2*(-l(i,1)+1),0,2*x(i,1),-1,0,0,0,0;
       0,2*l(i,1),2*x(i,2),0,-1,0,0,0;
       2*x(i,1),2*x(i,2),0,0,0,1,0,0;
       -1,0,0,0,0,0,1,0;
       0,-1,0,0,0,0,0,1;
       0,0,z(i,1),0,0,l(i,1),0,0;
       0,0,0,z(i,2),0,0,l(i,2),0;
       0,0,0,0,z(i,3),0,0,l(i,3)];
    p=-J\r(i,:)';
    z(i+1,:)=z(i,:)+a(i)*p(6:8)';
    l(i+1,:)=l(i,:)+a(i)*p(3:5)';
    while min(z(i+1,:))<=0
        a(i)=a(i)/2;
        z(i+1,:)=z(i,:)+a(i)*p(6:8)';
        l(i+1,:)=l(i,:)+a(i)*p(3:5)';
    end
    x(i+1,:)=x(i,:)+a(i)*p(1:2)';
    f(i+1)=2*(x(i+1,2))-x(i+1,1)^2-mu*e*log(z(i+1,1:3))'-l(i+1,:)*[1-x(i+1,1)^2-x(i+1,2)^2-z(i+1,1),x(i+1,1)-z(i+1,2),x(i+1,2)-z(i+1,3)]';
    while 1-x(i+1,1)^2-x(i+1,2)^2<0 || x(i+1,1)<0 || x(i+1,2)<0 || f(i+1)>f(i)-0.1*p'*r(i,:)'*a(i)
        a(i)=a(i)/2;
        x(i+1,:)=x(i,:)+a(i)*p(1:2)';
        z(i+1,:)=z(i,:)+a(i)*p(6:8)';
        l(i+1,:)=l(i,:)+a(i)*p(3:5)';
        f(i+1)=2*(x(i+1,2))-x(i+1,1)^2-mu*e*log(z(i+1,1:3))'-l(i+1,:)*[1-x(i+1,1)^2-x(i+1,2)^2-z(i+1,1),x(i+1,1)-z(i+1,2),x(i+1,2)-z(i+1,3)]';
    end
    if mu>tol
        mu=mu/2;
    end
    if (norm(x(i+1,:)-x(i,:))<tol*(1+norm(x(i+1,:))) && norm(r(i,:))<tol)
        break;
    end
end