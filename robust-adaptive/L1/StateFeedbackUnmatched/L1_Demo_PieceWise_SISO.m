% Demo of L1 adaptive control
% Unmatched uncertaintes with nonlinear dynamics
% Cut down SISO version
% Dr Hovakimyan Text Section 3.2.4
% Using aspects of modified control law from Banerjee et. al. 2016
% Joe Chai June 2019

%%
clearvars
clear all
n = 2;
m = 1;
% define system
% Am = [-10 -40;
%        1 0]; 
% Bm = [1;
%       0]; 
% c = [0 1];
Am = [0 1;-1 -1.4]; Bm = [0;1]; 
c = [1 0];
kg = -inv(c*inv(Am)*Bm);
Bum = [1; 
       0];
B = [Bm Bum];
   
% Uncertainties
Ad = [0.4 .5;
      0.5 0.2];
omega = 0.8;
omega_u = 5;
lambda = 0.3;

% solve lyapunov equation
Q = eye(n);
P = lyap(Am',Q);
%Pbar = [P(1,2); P(2,2)];

% sim params
dt = 0.0005;
t = 0:dt:30;
x_pred = zeros(n,length(t));
x           = zeros(n,length(t));
x_best      = zeros(n,length(t)); 
x_adapt     = zeros(n,length(t)); 
x_unmatched = zeros(n,length(t));
x_no_adapt  = zeros(n,length(t));

u_raw       = zeros(m,length(t));
u           = zeros(m,length(t));
u_best      = zeros(m,length(t));
u_no_adapt  = zeros(m,length(t));

sigma1_est  = zeros(m,length(t));
sigma2_est = zeros(n-m,length(t));

% adaptive parameters
sigma1_max = 5; 
sigma2_max = 5;

% adaptive law matrices
dt_adapt = dt;
PHI_adapt = inv(Am)*(expm(Am*dt_adapt) - eye(n));
eAmTs =  expm(Am*dt_adapt);


L = 10; % |5| + |5| (1 norm of vector)
%%
% Derive Hu and Hm
sym s;
s = tf('s');
Hu = c*inv(s*eye(n) - Am)*Bum; % unmatched transmission
num = Hu.Numerator{1};
den = Hu.Denominator{1};
[A_Hu,B_Hu,C_Hu,D_Hu] = tf2ss(num,den);

Hm_i = (1/-(c*inv(Am)*Bm)); % inverse of matched transmission

% Low pass filters C1 and C2
wk = 30;

%% run simulation
for i=3:length(t)
    r = sin(2*t(i)/pi);

    % adaptive law
    xtilde = x_pred(:,i-1) - x(:,i-1);
    sigma_est = -eye(2)*inv(B)*inv(PHI_adapt)*eAmTs*xtilde;
    sigma1_est(i-1) = sigma_est(1);
    sigma2_est(i-1) = sigma_est(2);
    
    % control law
    eta_1 = sigma1_est(i-1);
    eta_2 = sigma2_est(i-1);
    eta_2m = Hm_i * C_Hu * x_unmatched(:,i-1);
    
    u_raw(i-1) = (eta_1 + eta_2m - kg*r);
    u(i-1) = LPF_FO(-u_raw(i-1),u(i-2),dt,wk);
    u_best(:,i-1) = kg*r;
    
    % compute relevant dynamics
    xdot_pred      = Am*x_pred(:,i-1) ...
                     + Bm*(u(i-1) + sigma1_est(i-1)) ...
                     + Bum*(sigma2_est(i-1));
    x_best_dot     = Am*x_best(:,i-1) + Bm*kg*r;
    xdot          = (Am + Ad)*x(:,i-1) + Bm*omega*u(:,i-1) + [0.5*tanh(x(2,i-1))*x(2,i-1); 0.4*sin(3*x(1,i-1))*exp(x(2,i-1))];
    xdot_unmatched    = A_Hu*x_unmatched(:,i-1) + B_Hu*eta_2;
    xdot_no_adapt  = (Am + Ad)*x_no_adapt(:,i-1) + Bm*omega*kg*r + [0.5*tanh(x(2,i-1))*x(2,i-1); 0.4*sin(3*x(1,i-1))*exp(x(2,i-1))];
    
    % integrating/updating 
    x_pred(:,i)      = x_pred(:,i-1) + xdot_pred*dt;
    x_best(:,i)      = x_best(:,i-1) + x_best_dot*dt;
    x(:,i)           = x(:,i-1) + xdot*dt;
    x_unmatched(:,i) = x_unmatched(:,i-1) + xdot_unmatched*dt;
    x_no_adapt(:,i)  = x_no_adapt(:,i-1) + xdot_no_adapt*dt;
    
end
%%
figure('visible','on');
%subplot(2,1,1);
%plot(t,x_pred(1,:),'b--');
plot(t, x_best(1,:), 'r:');
hold on
plot(t, sin(2*t/pi), 'k-');
plot(t, x(1,:), 'm--');
plot(t, x_no_adapt(1,:), 'g--');
legend('y1 best','r1', 'y1 adapt','y1 no adapt');

figure('visible','on');
plot(t,sigma1_est);
hold on
plot(t,sigma2_est);
legend('sigma1 est','sigma2 est');

figure('visible','on');
plot(t,u_best(1,:));
hold on 
plot(t,u(1,:));
legend('u1 best','u');

figure('visible','on');
plot(t,x_pred - x);
legend('tracking error');

%% Analysis
% I = eye(2);
% s = tf('s');
% Hs = inv(s*I - Am)*Bm;
% Cs = wk/(s+wk);
% Gm = Hm*(1-Cs);
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

