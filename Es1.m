z1=2;
z2=2;
z3=2;
l1=2;
l2=2;
l3=2;
x1=1;
x2=1;
mu=0.2;
alpha=0.7;
eta=1;
gamma=1;
r=zeros(8,1);
while norm(r)<10e-6 
r1=2*(x1-4)-l1+l2;
r2=2*x2-l1+l3;
r3=2-x1-x2-z1;
r4=x1-z2;
r5=x2-z3;
r6=z1*l1-mu;
r7=z2*l2-mu;
r8=z3*l3-mu;
r=[r1,r2,r3,r4,r5,r6,r7,r8]';
J=[         2 , 0 , -1, 1 , 0 , 0 , 0 , 0;
            0 , 2 , -1, 0 , 1 , 0 , 0 , 0;
            1 , 1, 0 , 0 , 0 , -1, 0 ,0 ;
            1 , 0 , 0 , 0 , 0 , 0 , -1, 0;
            0 , 1 , 0 , 0 , 0 , 0 , 0 ,-1;
            0 , 0 , z1, 0 , 0 , l1, 0 , 0;
            0 , 0 , 0 , z2, 0 , 0 , l2, 0;
            0 , 0 , 0 , 0 , z3, 0 , 0 ,l3;
];
delta=J\r;
x1=x1+alpha*delta(1);
x2=x2+alpha*delta(2);
l1=l1+eta*delta(3);
l2=l2+eta*delta(4);
l3=l3+eta*delta(5);
z1=z1+gamma*delta(6);
z2=z2+gamma*delta(7);
z3=z3+gamma*delta(8);


end
