clear all
clc

%Parameters

M = 1.5;    % Mass of robot 
m = 0.046;  % Mass of wheels
b = 0.1;  % Viscous friction
MOI = 0.00414;% MOI of pendulum
g = 9.81;  %gravity
COG = 0.125;  %COG

%State-Space Matrices

den = MOI*(M+m)+M*m*COG^2; %denominator for the A and B matrices

A = [0      1              0           0;
     0 -(MOI+m*COG^2)*b/den  (m^2*g*COG^2)/den   0;
     0      0              0           1;
     0 -(m*COG*b)/den       m*g*COG*(M+m)/den  0]      %A-matrix
 
 
B = [     0;
     (MOI+m*COG^2)/den;
          0;
        m*COG/den]           %B-matrix
    
C = [1 0 0 0;
     0 0 1 0]
D = [0;
     0]                    %C-matrix

%Weighting Matrices

 Q = [1000     0     0     0
     0     0     0     0
     0     0     1000     0
     0     0     0     0]           %Q, weighting matrix
 
R = 1                               %R, weighting matrix

%gain calculation

Kr = lqr(A,B,Q,R)           %gain matrix