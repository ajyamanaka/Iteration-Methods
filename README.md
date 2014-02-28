Iteration-Methods
=================
function [ Iteration ] = myiterate( A,b,n,x0,w,s,tol )
%This function will take a matrix A of size nxn and vector b of size nx1
%and solve for the Jacobi iteration, Gauss Seidell method, and SOR method.
%x0 is the initial iterative vector. w is the omega factor in SOR.The 
%function will give a matrix where the columns of the matrix represent each 
%iteration vector.If s=1, the Jacobi iteration will be given. 
%If s=2,the Gauss Seidel will be given. 
%If s=3, SOR. For any other s, the function will output a 0. The tol is the
%tolerance of the solution given. The solution will give you an nx2 matrix.
%The first column will have the number of iterations and the second will
%have the iteration vector solution.
%Author Alan Yamanaka 

C=A\b;

%Jacobi

D=zeros(n,n);
for i=1:n
    D(i,i)=A(i,i);
end
R=A-D;
x=x0;
p=0;
while abs(C-x)>tol
    x=D^-1*(b-R*x);
    p=p+1;
end
P=zeros(n,1);
P(1,1)=p;



%Gauss Seidel
L=tril(A);
U=triu(A,1);
y=x0;
q=0;
while abs(C-y)>tol
    y=L^-1*(b-U*y);
    q=q+1;
end
Q=zeros(n,1);
Q(1,1)=q;

%SOR
L1=tril(A,-1);
U1=triu(A,1);
z=x0;
t=0;
while abs(C-z)>tol
    z=(D+w*L1)^-1*(w*b-(w*U1+(w-1)*D)*z);
    t=t+1;
end
T=zeros(n,1);
T(1,1)=t;


if s==1
    Iteration=[P x];
elseif s==2
        Iteration=[Q y];
elseif s==3
    Iteration=[T z];
    else Iteration=0;
end

        
       
    


end



Jacobi Gauss Seidel and SOR Methods Code
