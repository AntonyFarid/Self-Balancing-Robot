clear all
clc

% parameters of the SBR

M = 1.5;              % mass of the chassis, rods and shelves
m1 = 0.023;         % mass of wheel  
m= m1*2;            %mass of both wheels
b = 0.1;            % estimate of viscous friction coefficient (N-m-s)
g = 9.81;            %acceleration due to gravity (m/s^2)
l = 0.125;          %length to pendulum center of mass
r = 0.05;           %radius of wheel
d = 0.3;            %distance between wheels
Ir = (m1*r^2)/3;        % radial MOI of wheels
Ia = Ir*2;    % axial MOI of wheels
Ix   = 0.007402;    %x-axis MOI of pendulum      
Iy   = 0.00414;  %y-axis MOI of pendulum    
Iz   = 0.00234;    %z-axis MOI of pendulum  
kt = 0.658;          %torque constant 
ke = 1;          %back emf constant
R= 1;            %internal resistance of motor


%state-space eaquation

a_42_num = -(M^2)*(l^2)*g;                                          %numerator of the a42 component of the A matrix
a_52_num = ((M^2)*l*g)+(2*M*l*g*(m+(Ia/(r^2))));                    %numerator of the a52 component of the A matrix
a_44_num = 2*kt*ke*(((M*l^2)/(r^2))+(M*l/r)+(Iy/(r^2)));            %numerator of the a44 component of the A matrix
a_54_num = 2*kt*ke*(((2*Ia)/(r^2))+(M+2*m)+((l*M)/r));              %numerator of the a54 component of the A matrix
den= (M*Iy)+2*(Iy+M*(l^2))*(m+(Ia/(r^2)));                          %denomentator shared by all equations
a_42 = a_42_num/den;                                        % equation for a42
a_52 =a_52_num/den;                                         % equation for a52
a_44 = a_44_num/(R*den);                                    % equation for a44
a_54 = a_54_num/(R*r*den);                                  % equation for a54

A = [0  0       0   1       0   0;
    0   0       0   0       1   0;
    0   0       0   0       0   1;
    0   a_42    0   a_44    0   0;
    0   a_52    0   a_54    0   0;
    0   0       0   0       0   0]              %A-Matrix



b_41_num  =  (kt/R)*(((Iy+(M*l^2))/r)+M*l);
b_51_num  =(kt/R)*(((-M*(r+l))/r)-(2*(m+(Ia/(r^2))))) ;
b_61  = (kt*(d/(2*r)))/(R*(Iz+((d^2)/(2*r))*((m*r)+(Ia/r)))) ;


b_41 = b_41_num/den;
b_51 = b_51_num/den;
b_42 = b_41;
b_52 = b_51;
b_62 = -b_61;


B = [0  0;
    0   0;
    0   0;
    b_41  b_42;
    b_51  b_52;
    b_61  b_62]                     %B-Matrix

C   = eye(6)                        %C-Matrix
D   = zeros(6,2)                    %D-Matrix

SBR = ss(A,B,C,D);                  %State-Space model of the system


%LQR weighting matrices

Q   = diag([1000 5000 100 1000 1 0])           %Q-matrix
 
R   = diag([0.01 0.01])                         %R-matrix



%Gain matrix


K = lqr(A,B,Q,R)                                %K-gain Matrix

fb= feedback(SBR,K);
step(fb,'r')
