
% ------------------------------
% THIRD EXERCISE 
% WRITTEN BY MEHRAB ZAMANIAN
% ------------------------------
clear
clc
syms  E11 E22 G12 V12 t1 t2 theta t
% calculate stiffness matrix in y axis
ply=4;
thickness=[t/2-t1 t1 t1 t/2-t1];
V21=E22*V12/E11;
q=Q(E11,E22,G12,V12,V21);
z(1)=t1;z(2)=0;z(3)=t1;z(4)=t/2;
the(1)=theta;the(2)=-theta;the(3)=-theta;the(4)=theta;
for i=1:ply
T(:,:,i)=trans(the(i));
qbar(:,:,i)=T(:,:,i)*q*T(:,:,i)';
A(:,:,i)=qbar(:,:,i)*thickness(i);
if i==1
D(:,:,i)=1/3*qbar(:,:,i)*(z(i)^3-(-t/2)^3);
B(:,:,i)=1/2*qbar(:,:,i)*(z(i)^2-(-t/2)^2);
else
   D(:,:,i)=1/3*qbar(:,:,i)*(z(i)^3-z(i-1)^3);
   B(:,:,i)=1/2*qbar(:,:,i)*(z(i)^2-z(i-1)^2);
end
end
A_global=A(:,:,1)+A(:,:,2)+A(:,:,3)+A(:,:,4);
B_global=B(:,:,1)+B(:,:,2)+B(:,:,3)+B(:,:,4);
D_global=D(:,:,1)+D(:,:,2)+D(:,:,3)+D(:,:,4);
stiff(1:3,1:3)=A_global;
stiff(4:6,4:6)=D_global;
stiff(1:3,4:6)=B_global;
stiff(4:6,1:3)=B_global;
% solve equation D16=0 to obtain t1
    eqn=D_global(1,3)==0;






%DEFINE STIFFNESS MATRIX FUNCTION
function q=Q(E11,E22,G12,V12,V21)
q(1,1)=E11/(1-V21*V12);q(2,2)=E22/(1-V21*V12);q(3,3)=G12...
    ;q(2,1)=V21*E11/(1-V21*V12);q(1,2)=V12*E22/(1-V21*V12);
end
%DEFINE TRANSFORMATION MATRIX FUNCTION
function T=trans(theta)
M=cos(theta);
N=sin(theta);
T=[M^2 N^2 -2*M*N;N^2 M^2 2*M*N;N*M -M*N M^2-N^2];
end