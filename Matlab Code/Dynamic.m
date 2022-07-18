clc;
clear all;
close all;
%our data
w2 = -1;
a2 = 0;

r2 = 0.18; r1 = 0.5; h = 0.38;

% our unknown variables are r3, r4, p , teta4
%first loop
    %r2 * cos(teta2) - r3 * cos(teta4) = 0;
    %r2 * sin(teta2) - r3 * sin(teta4) + r1 = 0;
%second loop 
    %p - r4 * cos(teta4) = 0;
    %r1 + h - r4 * sin(teta4) = 0;

%we concider (r3, r4, p , teta4) as (x(1),x(2),x(3),x(4))
%the initial value of teta2
teta2 = pi/2;
%the initial values are achieved by solving equation with (teta2 = pi/2 rad)
firstGuess = [0.68;0.88;0; pi/2];
t = linspace(0,2*pi,361);
ans1 = [0;0;0;0];    
i = 1;
while teta2 >= -3*pi/2
    F1 = @(x)[ r2 * cos(teta2) - x(1) * cos(x(4));
               r2 * sin(teta2) - x(1) * sin(x(4)) + r1 ;
               x(3) - x(2) * cos(x(4));
               r1 + h - x(2) * sin(x(4)) ];
     x = fsolve(F1,firstGuess);
     ans1(:,i) = x;
     firstGuess = x;
     i = i + 1;
     teta2 = teta2 - 2*pi/360;
end

%solve the equations of velocity and angular velocity
%unknown variables are  r3_d , r4_d ,p_d , teta4_d ( v(1) v(2) v(3) v(4))
teta2 = pi/2;
ans2 = [0;0;0;0];    
k = 1;
while teta2 >= -3*pi/2
    V_coefficient = [-cos(ans1(4,k)), 0, 0, ans1(1,k) * sin(ans1(4,k));
                -sin(ans1(4,k)), 0, 0, -ans1(1,k) *cos(ans1(4,k));
                0, -cos(ans1(4,k)), 1, ans1(2,k)*sin(ans1(4,k));
                0, -sin(ans1(4,k)), 0, -ans1(2,k)*cos(ans1(4,k))];

    V_answer = [r2 * w2 *sin(teta2);
                -r2 * w2 *cos(teta2);
                 0;
                 0];
    ans2(:,k) = linsolve(V_coefficient,V_answer);
    k = k + 1;
    teta2 = teta2 - 2*pi/360;
          
end

%solve the equations of acceleration and angular accelration)
%unknown variables are  r3_dd , r4_dd ,p_dd , teta4_dd
teta2 = pi/2;
ans3 = [0;0;0;0];    
k = 1;
while teta2 >= -3*pi/2
    A_coefficient = [-cos(ans1(4,k)),0,0, ans1(1,k)*sin(ans1(4,k));
                    -sin(ans1(4,k)), 0, 0, -ans1(1,k)*cos(ans1(4,k));
                    0, -cos(ans1(4,k)), 1, ans1(2,k)*sin(ans1(4,k));
                    0, -sin(ans1(4,k)), 0, -ans1(2,k)*cos(ans1(4,k))];
	
	A_answer = [r2*a2*sin(teta2) + r2 *(w2)^2 *cos(teta2) - 2*ans2(1,k)*ans2(4,k)*sin(ans1(4,k)) - ans1(1,k)*(ans2(4,k)^2) * cos(ans1(4,k));
                -r2*(a2)*cos(teta2) + r2 *(w2)^2 *sin(teta2) + 2*ans2(1,k)*ans2(4,k)*cos(ans1(4,k)) - ans1(1,k)*(ans2(4,k)^2) * sin(ans1(4,k));
                -2*ans2(2,k)*ans2(4,k)*sin(ans1(4,k)) - ans1(2,k)*(ans2(4,k)^2)*cos(ans1(4,k));
                2*ans2(2,k)*ans2(4,k)*cos(ans1(4,k)) - ans1(2,k)*(ans2(4,k)^2)*sin(ans1(4,k))]; 
    ans3(:,k) = linsolve(A_coefficient,A_answer);
    k = k + 1;
    teta2 = teta2 - 2*pi/360;
end

% plot r3 (velocity and acceleration and displacement)
figure (1)
subplot(3,1,1)
plot(t,ans3(1,:))
legend('acceleration');
xlabel('t')
ylabel('r3_dd(m/s^2)')
grid on

subplot(3,1,2)
plot(t,ans2(1,:))
legend('velocity');
xlabel('t')
ylabel('r3_d(m/s)')
grid on

subplot(3,1,3)
plot(t,ans1(1,:))
legend('displacement');
xlabel('t')
ylabel('r3(m)')
grid on

% plot r4 (velocity and acceleration and displacement)
figure (2)
subplot(3,1,1)
plot(t,ans3(2,:))
legend('acceleration');
xlabel('t')
ylabel('r4_dd(m/s^2)')
grid on

subplot(3,1,2)
plot(t,ans2(2,:))
legend('velocity');
xlabel('t')
ylabel('r4_d(m/s)')
grid on

subplot(3,1,3)
plot(t,ans1(2,:))
legend('displacement');
xlabel('t')
ylabel('r4(m)')
grid on

% plot p (velocity and acceleration and displacement)
figure (3)
subplot(3,1,1)
plot(t,ans3(3,:))
legend('acceleration');
xlabel('t')
ylabel('p_dd(m/s^2)')
grid on

subplot(3,1,2)
plot(t,ans2(3,:))
legend('velocity');
xlabel('t')
ylabel('p_d(m/s)')
grid on

subplot(3,1,3)
plot(t,ans1(3,:))
legend('displacement');
xlabel('t')
ylabel('p(m)')
grid on

% plot teta4 (velocity and acceleration )
figure (4)
subplot(3,1,1)
plot(t,ans3(4,:))
legend('acceleration');
xlabel('t')
ylabel('teta4_dd(rad/s^2)')
grid on

subplot(3,1,2)
plot(t,ans2(4,:))
legend('velocity');
xlabel('t')
ylabel('teta4_d(rad/s)')
grid on

subplot(3,1,3)
plot(t,ans1(4,:))
legend('displacement');
xlabel('t')
ylabel('teta4(rad)')
grid on



%%%%%%%%%  PART3  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now we want to do the force analysis
g = 9.81;
p = 7801;
r4max = 0.9432;
k = 1;
teta2 = pi/2;
forceAnalysis = zeros(11,1);
while teta2 >= -3*pi/2
    aG2x = -r2/2 * (w2)^2 * cos(teta2);
    aG2y = -r2/2 * (w2)^2 * sin(teta2);
    aG3x = -r2 * (w2)^2 * cos(teta2);
    aG3y = -r2 * (w2)^2 * sin(teta2);
    aG4x = -r4max/2 * ((ans2(4,k))^2 *cos(ans1(4,k)) + ans3(4,k)*sin(ans1(4,k)));
    aG4y = r4max/2 * (-(ans2(4,k))^2 * sin(ans1(4,k))+ ans3(4,k)*cos(ans1(4,k)));
    aG5x = ans3(3,k);
    aG5y = 0;
    I2 = 0.0034;
    I3 = 0.000583;
    I4 = 0.4908;
    m2 = 1.26;
    m4 = 6.62;
    m3 = 0.7;
    m5 = 0.7; 
%%%%%%%%%%%%%%%%%%%%%%  FAx  FAy   F12x  F12y   M2   F34  F14x  F14y  FB  F15 T34
    CoeficientMatrix = [-1   ,0    ,1    ,0     ,0   ,0   ,0   ,0    ,0   ,0  ,0;
                        0    ,1    ,0    ,1     ,0   ,0   ,0   ,0    ,0   ,0  ,0;
                        -r2/2 * sin(2*pi-teta2),r2/2 * cos(2*pi-teta2)  ,-r2/2 * sin(2*pi-teta2),-r2/2 * cos(2*pi-teta2),1,0,0,0,0,0,0;
                        1    ,0    ,0    ,0     ,0   ,-sin(ans1(4,k)),0 ,0 ,0 ,0 ,0;
                        0    ,-1   ,0    ,0     ,0   ,cos(ans1(4,k)) ,0 ,0 ,0 ,0 ,0;
                        0    ,0    ,0    ,0     ,0   ,0   ,0   ,0    ,0   ,0  ,1;
                        0    ,0    ,0    ,0     ,0   ,sin(ans1(4,k)),1,0,sin(ans1(4,k)),0,0;
                        0    ,0    ,0    ,0     ,0   ,-cos(ans1(4,k)),0,1,-cos(ans1(4,k)),0,0;
                        0    ,0    ,0    ,0     ,0   ,r4max/2 - ans1(1,k),r4max/2 *sin(ans1(4,k)),-r4max/2 *cos(ans1(4,k)) ,-r4max/2,0,-1;
                        0    ,0    ,0    ,0     ,0   ,0,0,0,-sin(ans1(4,k)),0,0;
                        0    ,0    ,0    ,0     ,0   ,0,0,0,cos(ans1(4,k)),1,0];
                    
     AnswerMatrix = [(m2)*(aG2x);
                     (m2)*g + (m2)*(aG2y);
                     (I2)*(a2);
                     (m3)*(aG3x);
                     (m3)*g + (m3)*(aG3y);
                     (I3)*ans3(4,k);
                     (m4)*(aG4x);
                     (m4)*g + (m4)*(aG4y);
                     (I4)*ans3(4,k);
                     (m5)*(aG5x);
                     (m5)*g + (m5)*(aG5y)];
                 
    ANS = linsolve(CoeficientMatrix,AnswerMatrix);
    forceAnalysis(:,k) = ANS;
    k = k + 1;
    teta2 = teta2 - 2*pi/360;          
end
% plot FAx
figure (5)
plot(t,forceAnalysis(1,:))
title('FAx')
xlabel('t')
ylabel('N')
grid on

% plot FAy
figure (6)
plot(t,forceAnalysis(2,:))
title('FAy')
xlabel('t')
ylabel('N')
grid on

% plot F12x
figure (7)
plot(t,forceAnalysis(3,:))
title('F12x')
xlabel('t')
ylabel('N')
grid on

% plot F12y
figure (8)
plot(t,forceAnalysis(4,:))
title('F12y')
xlabel('t')
ylabel('N')
grid on

% plot M2
figure (9)
plot(t,forceAnalysis(5,:))
title('M2')
xlabel('t')
ylabel('N.m')
grid on

% plot F34
figure (10)
plot(t,forceAnalysis(6,:))
title('F34')
xlabel('t')
ylabel('N')
grid on

% plot F14x
figure (11)
plot(t,forceAnalysis(7,:))
title('F14x')
xlabel('t')
ylabel('N')
grid on

% plot F14y
figure (12)
plot(t,forceAnalysis(8,:))
title('F14y')
xlabel('t')
ylabel('N')
grid on

% plot FB
figure (13)
plot(t,-forceAnalysis(9,:))
title('FB')
xlabel('t')
ylabel('N')
grid on

% plot F15y
figure (14)
plot(t,forceAnalysis(10,:))
title('F15y')
xlabel('t')
ylabel('N')
grid on

% plot T34
figure (15)
plot(t,forceAnalysis(11,:))
title('T34')
xlabel('t')
ylabel('N.m')
grid on