function y= Conjugate_gradient(A,F,x0,ep)
%Conjugate_gradient method for solving system of linear equations
r0=F-A*x0;
p0=r0;
alpha0=dot(r0,p0)/dot((A*p0),p0);
x1=x0+alpha0*p0;
r1=F-A*x1;
n=1;
while norm(x1-x0,inf)>=ep
    beta0=-dot(r1,A*p0)/dot(p0,A*p0);
    p1=r1+beta0*p0;
    alpha1=dot(r1,p1)/dot((A*p1),p1);
    x2=x1+alpha1*p1;
    r2=F-A*x2;
    r1=r2;
    p0=p1;
    X=x2;
    x0=x1;
    x1=x2;
    n=n+1;
end
y=X;
disp('Number of Iterations');
disp(n);
end