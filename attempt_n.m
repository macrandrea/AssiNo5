maxit=1000;
tol=1e-8;
r=ones(maxit,8);
x=ones(maxit,2);
l=ones(maxit,3);
z=ones(maxit+1,3);
a=ones(maxit,1);
f=ones(maxit+1,1);
mu=1;
x(1,1:2)=rand(1,2)*2;
if sum(x(1,:))>2
    x(1,:)=x(1,:)-(sum(x(1,:))-2);
end
l(1,1:3)=rand(1,3);
z(1,1:3)=[2-x(1,1)-x(1,2),x(1,1),x(1,2)];
e=[1,1,1];
f(1)=(x(1,1)-4)^2-x(1,2)^2-mu*e*log(z(1,1:3))'+l(1,:)*[2-x(1,1)-x(1,2)-z(1,1),x(1,1)-z(1,2),x(1,2)-z(1,3)]';
for i=1:maxit
    r(i,1)=2*(x(i,1)-4)-l(i,1)+l(i,2);%deriv Lagr di x1    
    r(i,2)=2*x(i,2)-l(i,1)+l(i,3);    %deriv Lagr di x2     
    r(i,3)=2-x(i,1)-x(i,2)-z(i,1);    %deriv Lagr di l1  
    r(i,4)=x(i,1)-z(i,2);         %deriv Lagr di l2 
    r(i,5)=x(i,2)-z(i,3);         %deriv Lagr di l3 
    r(i,6)=-l(i,1)-mu/z(i,1);      %deriv Lagr di z1  
    r(i,7)=-l(i,2)-mu/z(i,2);      %deriv Lagr di z2     
    r(i,8)=-l(i,3)-mu/z(i,3);      %deriv Lagr di z3
    J=[2,0,-1,1,0,0,0,0;
       0,2,-1,0,1,0,0,0;
       -1,-1,0,0,0,-1,0,0;
       1,0,0,0,0,0,-1,0;
       0,1,0,0,0,0,0,-1;
       0,0,-1,0,0,mu/z(i,1)^2,0,0;
       0,0,0,-1,0,0,mu/z(i,2)^2,0;
       0,0,0,0,-1,0,0,mu/z(i,3)^2];
    p=-J\r(i,:)';
    z(i+1,:)=z(i,:)+a(i)*p(6:8)';
    l(i+1,:)=l(i,:)+a(i)*p(3:5)';
    while min(z(i+1,:))<=0
        a(i)=a(i)/2;
        z(i+1,:)=z(i,:)+a(i)*p(6:8)';
        l(i+1,:)=l(i,:)+a(i)*p(3:5)';
    end
    x(i+1,:)=x(i,:)+a(i)*p(1:2)';
    f(i+1)=(x(i+1,1)-4)^2-(x(i+1,2))^2-mu*e*log(z(i+1,:))'+l(i+1,:)*[2-x(i+1,1)-x(i+1,2)-z(i+1,1),x(i+1,1)-z(i+1,2),x(i+1,2)-z(i+1,3)]';
    while x(i+1,1)<0 || x(i+1,2)<0 || x(i+1,1)+x(i+1,2)>2 || f(i+1)>f(i)
        a(i)=a(i)/2;
        x(i+1,:)=x(i,:)+a(i)*p(1:2)';
        z(i+1,:)=z(i,:)+a(i)*p(6:8)';
        l(i+1,:)=l(i,:)+a(i)*p(3:5)';
        f(i+1)=(x(i+1,1)-4)^2-(x(i+1,2))^2-mu*e*log(z(i+1,:))'+l(i+1,:)*[2-x(i+1,1)-x(i+1,2)-z(i+1,1),x(i+1,1)-z(i+1,2),x(i+1,2)-z(i+1,3)]';
    end
    mu=mu/2;
    if (norm(x(i+1,:)-x(i,:))<tol*(1+norm(x(i+1,:))))
        break;
    end
end