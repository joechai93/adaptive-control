% Demo of L1 adaptive control
% Dr Hovakimyan Text Section 2.2.7
% Joe Chai June 2019

clearvars
clear all
% define system
Am = [0 1;-1 -1.4]; B = [0;1]; 
c = [1; 0];
kg = -1/(c'*inv(Am)*B);

% case 1
omega = 1; 
theta = @(t) [2+cos(pi*t); 2+0.3*sin(pi*t)+0.2*cos(2*t)];
sigma = @(t) sin(pi*t);

% solve lyapunov equation
Q = eye(2);
P = lyap(Am',Q);
%Pbar = [P(1,2); P(2,2)];

% sim params
dt = 0.00001;
t = 0:dt:25;
x_pred = zeros(2,length(t));
x      = zeros(2,length(t));
u      = zeros(length(t),1);
u_raw  = zeros(length(t),1);
theta_est  = zeros(2,length(t)); theta_est(:,1) = theta(0); theta_est(:,2) = theta(0);
sigma_est  = zeros(1,length(t)); sigma_est(1) = sigma(0);sigma_est(2) = sigma(0);
omega_est  = zeros(1,length(t)); omega_est(1) = 1;omega_est(2) = 1;

% adaptive parameters
gam = 50000;
theta_max = 5;
sigma_max = 50;
omega0 = [0.1 2];
omega_mid = (omega0(1) + omega0(2))/2;
omega_d_max = abs(omega0(2) - omega_mid);

L = 10; % |5| + |5| (1 norm of vector)
% Low pass filter C(s) 
% 0.5 hz * k=60
wk = 30;

tic
% run simulation
for i=3:length(t)
    r = cos(2*t(i)/pi);
    % control law
    eta_est = theta_est(:,i-1)'*x(:,i-1) + sigma_est(i-1);
    u_raw(i-1) = (eta_est - kg*r)/omega_est(i-1);
    u(i-1) = LPF_FO(-u_raw(i-1),u(i-2),dt,wk);
    % adaptive law
    xtilde = x_pred(:,i-1) - x(:,i-1);
    theta_estdot  = gam * Proj(theta_est(:,i-1), -xtilde'* P * B * x(:,i-1),theta_max);
    sigma_estdot  = gam * Proj(sigma_est(i-1), -xtilde'* P * B ,sigma_max);
    omega_estdot  = gam * Proj((omega_est(i-1)-omega_mid), -xtilde'* P * B * u(i-1),omega_d_max);
    % compute relevant dynamics
    xdot_pred = Am*x_pred(:,i-1) + B*(theta_est(:,i-1)'*x(:,i-1) + omega_est(i-1)*u(i-1) + sigma_est(i-1));
    xdot     = Am*x(:,i-1) + B*(theta(t(i-1))'*x(:,i-1) + omega*u(i-1) + sigma(t(i-1)));
    % integrating/updating 
    x_pred(:,i) = x_pred(:,i-1) + xdot_pred*dt;
    x(:,i)     = x(:,i-1) + xdot*dt;
    theta_est(:,i)  = theta_est(:,i-1) + theta_estdot*dt;
    sigma_est(i)  = sigma_est(i-1) + sigma_estdot*dt;
    omega_est(i)  = omega_est(i-1) + omega_estdot*dt;
    if i > length(t)/4 && i < length(t)/2
        disp('25 percent done');
    elseif i > length(t)/2 && i < 3*length(t)/4
        disp('50 percent done');
    elseif i > 3*length(t)/4 && i < length(t)
        disp('75 percent done');
    end
end
%disp('100 percent done');
toc


figure('visible','on');
plot(t,x_pred(1,:),'b--');
hold on
plot(t,x(1,:),'r--');
plot(t,cos(2*t/pi));
legend('predicted','true','ref');

figure('visible','on');
plot(t,theta_est(1,:));
hold on
plot(t,theta_est(2,:));
theta_plt = theta(t);
plot(t,theta_plt(1,:),'o');
plot(t,theta_plt(2,:),'*');
legend('theta1 est','thets2 est','theta1 true','theta2 true');

figure('visible','on');
plot(t,sigma_est(:));
hold on
plot(t,sigma(t),'o');
legend('sigma est','sigma true');

figure('visible','on');
plot(t,omega_est(:));
hold on
plot([0 t(end)],omega);
legend('omega est','omega true');

figure('visible','on');
plot(t,u);
hold on
plot(t,-u_raw,'k:');
legend('u lpf','-u');

%% Analysis
I = eye(2);
s = tf('s');
Hs = inv(s*I - Am)*B;
Cs = wk/(s+wk);
Gs = Hs*(1-Cs);

[G1_num, G1_den] = tfdata(Gs(1),'v');
[G2_num, G2_den] = tfdata(Gs(2),'v');
 
[A1,B1,C1,D1] = tf2ss(G1_num,G1_den); 
[A2,B2,C2,D2] = tf2ss(G2_num,G2_den);

G1_ss = ss(A1,B1,C1,D1);
G2_ss = ss(A2,B2,C2,D2);

% use G1_ss as c = [1 0]
disp('lambda is')
L*(l1norm(G1_ss)+l1norm(G2_ss))


