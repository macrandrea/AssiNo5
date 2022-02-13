clear
clc
maxit=100;
tol=1e-8;
mat=zeros(maxit+1,2);
alpha=ones(maxit,1);
f=zeros(maxit+1,1);
mat(1,:)=[8,0.8];
x=mat(1,1);
y=mat(1,2);
f(1)=(3/2-x*(1-y))^2+(9/4-x*(1-y^2))^2+(21/8-x*(1-y^3))^2;
for i=1:maxit
    H=2*[(1 - y^3)^2 + (1 - y^2)^2 + (1 - y)^2, -2* x*y* (1 - y^2) + 2*y* (9/4 - x*(1 - y^2)) - 3* x *(1 - y^3)* y^2 + 3* y^2* (21/8 - x* (1 - y^3)) + 3/2 - x *(1 - y) - x *(1 - y);
        -2* x* (1 - y^2)* y + 2* y *(9/4 - x *(1 - y^2)) - 3* x* (1 - y^3)* y^2 + 3* y^2* (21/8 - x *(1 - y^3)) + 3/2 - x* (1 - y) - x *(1 - y), 9* x^2 *y^4 + 4* x^2* y^2 + x^2 + 6 *x* y* (21/8 - x *(1 - y^3)) + 2* x* (9/4 - x *(1 - y^2))];
    J=[(y - 1)* (8*x*(y^5+y^4+2*y^3-y-3)+21*y^2+39*y+51); x*(8*x*(3*y^5+2*y^3-3*y^2-y-1)+63*y^2+36*y+12)]/4;
    p=H\J;
    slope=p'*J;
    if slope>0
        while(true)
            mat(i+1,:) = mat(i,:)-alpha(i)*p';
            x=mat(i+1,1);
            y=mat(i+1,2);
            f(i+1)=(3/2-x*(1-y))^2+(9/4-x*(1-y^2))^2+(21/8-x*(1-y^3))^2;
            if f(i+1)<f(i)-0.1*alpha(i)*slope
                break
            else
                alpha(i)=alpha(i)/2;
            end
        end
    end
    if (norm(mat(i+1,:)-mat(i,:))<tol*(1+norm(mat(i+1,:))) && norm(J)<tol)
        break;
    end   
end