% Demo of L1 adaptive control
% Dr Hovakimyan Text Section 2.1.5
% Joe Chai May 2019

clearvars
clear all
% define system
A = [0 1;-1 -1.4]; B = [0;1]; 
theta = [4; -4.5];
c = [1; 0];
km = [1.5; 1.5]; % textbook problem uses km = [0;0];
Am = A - B*km'; 
% solve lyapunov equation
Q = eye(2);
P = lyap(Am',Q);
%Pbar = [P(1,2); P(2,2)];

% sim params
dt = 0.00005; % must make this very small or signals will explode
t = 0:dt:40;
x_pred = zeros(2,length(t));
x      = zeros(2,length(t));
x_no_adapt = zeros(2,length(t));
x_ref = zeros(2,length(t));

u      = zeros(length(t),1);
u_raw  = zeros(length(t),1);
theta_est  = zeros(2,length(t)); theta_est(:,2) = [0;0];
control_dt = dt;

% adaptive parameters
gam = 10; 
kg = -1/(c'*inv(Am)*B);
theta_max = 10;
L = theta_max *2;
% Low pass filter C(s) 
wc = 5;

xdot_pred = 0;
xdot = 0;

eta_alc = 4.0; %4.0;
% setup dynamics
%f     = @(x,u) Am*x + B*(theta'*x + u);

% run simulation
for i=3:length(t)
    r = 25*sin(0.5*t(i));
    % control law
    u_raw(i-1) = theta_est(:,i-1)'*x(:,i-1) - kg*r;
    u(i-1) = LPF_FO(-u_raw(i-1),u(i-2),control_dt,wc);
    % adaptive law
    xtilde = x_pred(:,i-1) - x(:,i-1);

    %adaptive_law  = @(theta_est,u) gam * Proj(theta_est, -xtilde'* P * B * x(:,i-1),theta_max); % main one
    
    theta_estdot = gam * Proj(theta_est(:,i-1), -xtilde'* P * B * x(:,i-1) - eta_alc*theta_est(:,i-1) ,theta_max); % main one
    
    % compute relevant dynamics
    %f_pred = @(x_pred,u) Am*x_pred + B*(theta_est(:,i-1)'*x(:,i-1) + u);
    %xdot_pred = f_pred(x_pred(:,i-1),u(i-1)); 
    xdot_pred = Am*x_pred(:,i-1) + B*(theta_est(:,i-1)'*x(:,i-1) + u(i-1));
    %xdot     = f(x(:,i-1),u(i-1)); 
    xdot = A*x(:,i-1) + B*(theta'*x(:,i-1) + u(i-1) -km'*x(:,i-1));
    %xdot_no_adapt = Am*x_no_adapt(:,i-1) + B*(theta'*x_no_adapt(:,i-1) +kg*r);
    xdot_ref = A*x_ref(:,i-1) + B*(-km'*x_ref(:,i-1)+ kg*25*sin(0.5*t(i-1)));
    
    % integrating/updating 
    x_pred(:,i) = x_pred(:,i-1) + xdot_pred*dt;
    x(:,i)     = x(:,i-1) + xdot*dt;
    %x_no_adapt(:,i) = x_no_adapt(:,i-1) + xdot_no_adapt*dt;
    x_ref(:,i) = x_ref(:,i-1) + xdot_ref*dt;
    theta_est(:,i)  = theta_est(:,i-1) + theta_estdot*dt;
    
    %x_pred(:,i) = runge_kutta_four(x_pred(:,i-1),f_pred,u(i-1),dt);%x_pred(:,i-1) + xdot_pred*dt;
    %x(:,i)     = runge_kutta_four(x(:,i-1),f,u(i-1),dt);%x(:,i-1) + xdot*dt;
    %theta_est(:,i)  = runge_kutta_four(theta_est(:,i-1),adaptive_law,u(i-1),dt);%theta_est(:,i-1) + theta_estdot*dt;
    
    % percentage complete
    if i < 0.25*length(t)
        disp('0%');
    elseif i > 0.25*length(t) && i < 0.5*length(t)
        disp('25%');
    elseif i > 0.5*length(t) && i < 0.75*length(t)
        disp('50%');
    else
        disp('75%');
    end
end

figure('visible','on');
plot(t,x(1,:),'r--');
hold on;
plot(t,25*sin(0.5*t),'g--');
plot(t,x_ref(1,:),'m--');
legend('with adapt','com','ref');

figure('visible','on');
plot(t,theta_est(1,:));
hold on
plot(t,theta_est(2,:));
legend('theta1 est','thets2 est');

figure('visible','on');
plot(t,u);
hold on
plot(t,-u_raw,'k:');
legend('u lpf','-u');

%% Analysis
I = eye(2);
s = tf('s');
Hs = inv(s*I - Am)*B;
Cs = wc/(s+wc);
Gs = Hs*(1-Cs);

[G1_num, G1_den] = tfdata(Gs(1),'v');
[G2_num, G2_den] = tfdata(Gs(2),'v');
 
[A1,B1,C1,D1] = tf2ss(G1_num,G1_den); 
[A2,B2,C2,D2] = tf2ss(G2_num,G2_den);

G1_ss = ss(A1,B1,C1,D1);
G2_ss = ss(A2,B2,C2,D2);

% use G1_ss as c = [1 0]
%disp('lambda is')
%L*(l1norm(G1_ss)+l1norm(G2_ss))

%L*L1NormTwoState(Gs)