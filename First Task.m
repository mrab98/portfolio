
% ------------------------------
% FIRST EXERCISE 
% WRITTEN BY MEHRAB ZAMANIAN
% ------------------------------
clear
clc
syms theta
% INPUT MATERIAL PROPERTIES
E11=164;E22=8.30;G12=2.10;V12=0.34;
V21=E22*V12/E11;
%CALCULATE STIFFNESS MATRIX
q=Q(E11,E22,G12,V12,V21);
%CALCULATE TRANSFORMATION MATRIX
T=trans(theta);
%CALCULATE QBAR
qbar=T*q*T';
%plot qbar vs theta graph
y1=double(subs(qbar(1,1),theta,0:0.1:pi));
y2=double(subs(qbar(1,2),theta,0:0.1:pi));
y3=double(subs(qbar(2,1),theta,0:0.1:pi));
y4=double(subs(qbar(2,2),theta,0:0.1:pi));
y5=double(subs(qbar(3,3),theta,0:0.1:pi));

subplot(5,1,1)
plot(0:0.1:pi,y1);
xlabel('0 \leq Theta \leq \pi')
ylabel('Q11(Gpa)')

subplot(5,1,2)
plot(0:0.1:pi,y2);
xlabel('0 \leq Theta \leq \pi')
ylabel('Q12(Gpa)')

subplot(5,1,3)
plot(0:0.1:pi,y3);
xlabel('0 \leq Theta \leq \pi')
ylabel('Q21(Gpa)')

subplot(5,1,4)
plot(0:0.1:pi,y4);
xlabel('0 \leq Theta \leq \pi')
ylabel('Q22(Gpa)')

subplot(5,1,5)
plot(0:0.1:pi,y5);
xlabel('0 \leq Theta \leq \pi')
ylabel('Q33(Gpa)')

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
T=[M^2 N^2 2*M*N;N^2 M^2 -2*M*N;-N*M M*N M^2-N^2];
end
