
% ------------------------------
% SECIND EXERCISE 
% WRITTEN BY MEHRAB ZAMANIAN
% ------------------------------
clear
clc
 % INPUT MATERIAL PROPERTIES
E11=164;E22=8.30;G12=2.10;V12=0.34;
V21=E22*V12/E11;
t=0.01;
% number of laminas
ply=input('number of laminas');
% orientation of laminas pi/2
for i=1:ply
    theta(i)=input('orientation of laminas(bottom to top)');
end
%CALCULATE STIFFNESS MATRIX OF LAMINA IN LOCAL DIRECTION
q=Q(E11,E22,G12,V12,V21);
%CALCULATE Z VALUES
h=ply*t;
for i=1:ply
    if i==1
    z(i)=-h/2+t;
    else
     z(i)=z(i-1)+t;
    end
end

%CALCULATE STIFFNESS MATRIX OF LAMINAS IN GOLOBAL DIRECTION
%CALCULATE TRANSFORMATION MATRIX IN DIFFERNT DEGREE
for i=1:ply
T(:,:,i)=trans(theta(i));
qbar(:,:,i)=T(:,:,i)*q*T(:,:,i)';
A(:,:,i)=qbar(:,:,i)*t;
if i==1
D(:,:,i)=1/3*qbar(:,:,i)*(z(i)^3-(-h/2)^3);
B(:,:,i)=1/2*qbar(:,:,i)*(z(i)^2-(-h/2)^2);
else
   D(:,:,i)=1/3*qbar(:,:,i)*(z(i)^3-z(i-1)^3);
   B(:,:,i)=1/2*qbar(:,:,i)*(z(i)^2-z(i-1)^2);
end
end
%laminate stiffness with respect to the global x-y coordinate system
A_global=sum(A,3);
B_global=sum(B,3);
D_global=sum(D,3);
stiff(1:6,1:6)=0;
stiff(1:3,1:3)=A_global;
stiff(4:6,4:6)=D_global;
stiff(1:3,4:6)=B_global;
stiff(4:6,1:3)=B_global;


%DEFINE STIFFNESS MATRIX FUNCTION
function q=Q(E11,E22,G12,V12,V21)
q(1:3,1:3)=0;
q(1,1)=E11/(1-V21*V12);q(2,2)=E22/(1-V21*V12);q(3,3)=G12...
    ;q(2,1)=V21*E11/(1-V21*V12);q(1,2)=V12*E22/(1-V21*V12);
end
%DEFINE TRANSFORMATION MATRIX FUNCTION
function T=trans(theta)
M=cos(theta);
N=sin(theta);
T=[M^2 N^2 -2*M*N;N^2 M^2 2*M*N;N*M -M*N M^2-N^2];
end