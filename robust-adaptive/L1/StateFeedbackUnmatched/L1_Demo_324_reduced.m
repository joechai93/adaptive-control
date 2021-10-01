% Demo of L1 adaptive control
% Unmatched uncertaintes without unmodelled nonlinear dynamics
% Dr Hovakimyan Text Section 3.2.4
% Joe Chai June 2019

clearvars
clear all
n = 3;
m = 2;
% define system
Am = [-1 0 0;
       0 0 1;
       0 -1 -1.8]; 
Bm = [1 0;
      0 0;
      1 1]; 
c = [1 0 0;
     0 1 0];
kg = -inv(c*inv(Am)*Bm);
Bum = [0;1;0];
% scenario 1
Ad = [0.2 -0.2 -0.3;
      -0.2 -0.2 0.6;
      -0.1 0 -0.9];
omega = [0.6 -0.2;
         0.2 1.2];
omega_u = 5;
lambda = 0.3;

% solve lyapunov equation
Q = eye(3);
P = lyap(Am',Q);
%Pbar = [P(1,2); P(2,2)];

% sim params
dt = 0.0001;
t = 0:dt:30;
x_pred = zeros(n,length(t));
x           = zeros(n,length(t)); x(2,1) = 1;
x_best      = zeros(n,length(t)); x_best(2,1) = 1;
x_adapt     = zeros(n,length(t)); x_adapt(2,1) = 1;

u           = zeros(m,length(t));
u_best      = zeros(m,length(t));

theta1_est  = zeros(m,length(t));
sigma1_est  = zeros(m,length(t));
theta2_est = zeros(n-m,length(t));
sigma2_est = zeros(n-m,length(t));
omega11_est  = zeros(1,length(t));
omega22_est  = zeros(1,length(t));
omega12_est  = zeros(1,length(t));
omega21_est  = zeros(1,length(t));


% adaptive parameters
gam = 50000;
K = [8 0;0 8];
theta1_max = norm([40;40]);
theta2_max = 40;
sigma1_max = norm([5; 5]);
sigma2_max = 5;
omegaoff = 0.2;
omegadiag = [0.25 3];
omegadiag_mid = (omegadiag(1) + omegadiag(2))/2;
omegadiag_d_max = abs(omegadiag(2) - omegadiag_mid);
omega_max = norm([3 0.2;0.2 3]);

L = 10; % |5| + |5| (1 norm of vector)
% Low pass filter D(s) 
x1_filter = zeros(5,length(t)); 
eta1_filter = zeros(1,length(t));
x2_filter = zeros(5,length(t));
eta2_filter = zeros(1,length(t));

A_lpf =  [-167      -10190     -278000    -2800000 0;
           1           0           0           0    0;
           0           1           0           0    0;
           0           0           1           0    0;
           0           0           0           1    0];
B_lpf =[ 1;
        0;
        0;
        0;
        0];

C_lpf = [0           0           0     0   2800000];

tic
% run simulation
for i=2:length(t)
    r = [sin(pi*t(i)/3);
         0.2+0.8*cos(pi*t(i)/6)];
         
    % control law
    eta_1 = theta1_est(:,i-1)*x(:,i-1) + sigma1_est(i-1);
    eta_2 = theta2_est(:,i-1)*x(:,i-1) + sigma2_est(i-1);
    Hm = 
    
    u_raw(i-1) = eta_1 + eta_2m - kg*r;
    u_best(:,i-1) = %kg*r;%LPF_FO(-u_raw(i-1),u(i-2),dt,wk);
    % adaptive law
    xtilde = x_pred(:,i-1) - x(:,i-1);
    theta1_estdot  = gam * Proj(theta1_est(:,i-1), -(xtilde'* P * Bm)' * norm(x(:,i-1),'Inf'), theta1_max);
    theta2_estdot  = gam * Proj(theta2_est(:,i-1), -(xtilde'* P * Bum)' * norm(x(:,i-1),'Inf'), theta2_max); 
    sigma1_estdot  = gam * Proj(sigma1_est(:,i-1), -(xtilde'* P * Bm)' , sigma1_max);
    sigma2_estdot  = gam * Proj(sigma2_est(:,i-1), -(xtilde'* P * Bum)' , sigma2_max);
    %omega_est = [omega11_est(i-1) omega12_est(i-1);
    %         omega21_est(i-1) omega22_est(i-1)];
    %omega_estdot  = gam * Proj([omega_est(:,1); omega_est(:,2)], -(xtilde'* P * Bm)' * u(:,i-1)', omega_max);
    % compute relevant dynamics
    xdot_pred      = Am*x_pred(:,i-1) + Bm*(theta1_est(:,i-1)*norm(x(:,i-1),'Inf') + omega*u(:,i-1) + sigma1_est(:,i-1)) ...
                     + Bum*(theta2_est(:,i-1)*norm(x(:,i-1),'Inf') + sigma2_est(:,i-1));
    x_best_dot     = Am*x_best(:,i-1) + Bm*(u_best(:,i-1));
    x_dot          = (Am + Ad)*x(:,i-1) + Bm*omega*(u_best(:,i-1));
    x1_filter_dot  = A_lpf*x1_filter(:,i-1) + B_lpf*eta_est(1);
    x2_filter_dot  = A_lpf*x2_filter(:,i-1) + B_lpf*eta_est(2); 
    % integrating/updating 
    x_pred(:,i) = x_pred(:,i-1) + xdot_pred*dt;
    x_best(:,i)     = x_best(:,i-1) + x_best_dot*dt;
    x(:,i)          = x(:,i-1) + x_dot*dt;
    x1_filter(:,i)  = x1_filter(:,i-1) + x1_filter_dot*dt;
    x2_filter(:,i)  = x2_filter(:,i-1) + x2_filter_dot*dt;
    eta1_filter(i) = C_lpf*x1_filter(:,i);
    eta2_filter(i) = C_lpf*x2_filter(:,i);
    
    theta1_est(:,i)  = theta1_est(:,i-1) + theta1_estdot*dt;
    theta2_est(i)  = theta2_est(:,i-1) + theta2_estdot*dt;
    sigma1_est(:,i)  = sigma1_est(i-1) + sigma1_estdot*dt;
    sigma2_est(i)  = sigma2_est(i-1) + sigma2_estdot*dt;
    
    %omega_est = omega_est + [omega_estdot(1) omega_estdot(3); omega_estdot(2) omega_estdot(4)] * dt;
    %omega11_est(i) = omega_est(1,1);
    %omega12_est(i) = omega_est(1,2);
    %omega21_est(i) = omega_est(2,1);
    %omega22_est(i) = omega_est(2,2);
    
%     if i > length(t)/4 && i < length(t)/2
%         disp('25 percent done');
%     elseif i > length(t)/2 && i < 3*length(t)/4
%         disp('50 percent done');
%     elseif i > 3*length(t)/4 && i < length(t)
%         disp('75 percent done');
%     end
end
%disp('100 percent done');
toc


figure('visible','on');
subplot(2,1,1);
%plot(t,x_pred(1,:),'b--');
plot(t, x_best(1,:), 'r:');
hold on
plot(t, sin(pi*t/3), 'k-');
plot(t, x(1,:), 'r--');
legend('y1 best','r1', 'y1 no adapt');
subplot(2,1,2);
plot(t, x_best(2,:), 'b:');
hold on;
plot(t, 0.2+0.8*cos(pi*t/6), 'k-');
plot(t, x(2,:), 'b--');
legend('y2 best','r2','y2 no adapt');

% figure('visible','on');
% plot(t,theta_est(1,:));
% hold on
% plot(t,theta_est(2,:));
% theta_plt = theta(t);
% plot(t,theta_plt(1,:),'o');
% plot(t,theta_plt(2,:),'*');
% legend('theta1 est','thets2 est','theta1 true','theta2 true');
% 
% figure('visible','on');
% plot(t,sigma_est(:));
% hold on
% plot(t,sigma(t),'o');
% legend('sigma est','sigma true');
% 
% figure('visible','on');
% plot(t,omega_est(:));
% hold on
% plot([0 t(end)],omega);
% legend('omega est','omega true');

figure('visible','on');
subplot(2,1,1)
plot(t,u_best(1,:));
hold on
plot(t,eta1_filter);
legend('u1 best');
subplot(2,1,2);
plot(t,u_best(2,:));
hold on;
plot(t,eta2_filter);
legend('u2 best');

%% Analysis
% I = eye(2);
% s = tf('s');
% Hs = inv(s*I - Am)*B;
% Cs = wk/(s+wk);
% Gs = Hs*(1-Cs);
% 
% [G1_num, G1_den] = tfdata(Gs(1),'v');
% [G2_num, G2_den] = tfdata(Gs(2),'v');
%  
% [A1,B1,C1,D1] = tf2ss(G1_num,G1_den); 
% [A2,B2,C2,D2] = tf2ss(G2_num,G2_den);
% 
% G1_ss = ss(A1,B1,C1,D1);
% G2_ss = ss(A2,B2,C2,D2);
% 
% % use G1_ss as c = [1 0]
% disp('lambda is')
% L*(l1norm(G1_ss)+l1norm(G2_ss))
% 

