% Demo of indirect MIAC using least squares
% for first order siso system
% Dr Nhan Nguyen Text on MRAC section 6.3
% Joe Chai 2019
% a unknown, b known, theta unknown
% RLS and gradient least squares both implemented (comment out)
clearvars
close all
% define system
a = 1; b = 2; theta = 0.2;
am = -1; bm = 1;
Phi = @(x) x*x; 
Psi = @(x) [x; Phi(x)];
% sim params
dt = 0.01;
t = 0:dt:40;
x_ref = zeros(length(t),1);
x     = zeros(length(t),1);
u     = zeros(length(t),1);
a_est      = zeros(length(t),1); a_est(1) = a;
theta_est  = zeros(length(t),1); 
Omega_est  = zeros(2,length(t));
e_no_adapt = zeros(length(t),1);
x_no_adapt = zeros(length(t),1);
R = 10*eye(2); %RLS
r11 = zeros(length(t),1); r11(1) = R(1,1);
r12 = zeros(length(t),1);
r21 = zeros(length(t),1); r21(1) = R(2,1);
r22 = zeros(length(t),1);

% adaptive parameters
gam_a = 10;
gam_f = 10;
P = 1;
Gamma = [gam_a 0;0 gam_f]; % Gradient LS
% true gains
kx_true = (am - a)/b;
kr_true = bm/b;

% run smulation
r = 0;
xdot = 0;
for i=2:length(t)
    if t(i) > 0.5
        r = sin(t(i));
    end
    e = x_ref(i-1) - x(i-1);
    % control law
    kx = (am - a_est(i-1))/b;
    kr = bm/b;
    ubar = kx*x(i-1) + kr*r;
    xdot_des = a_est(i-1)*x(i-1) + b*ubar;
    epsilon = xdot_des - xdot;
    u(i-1) = kx*x(i-1) + kr*r - theta_est(i-1)*Phi(x(i-1)); 
    u_no_adapt = kx_true*x_no_adapt(i-1) + kr_true * r;
    % adaptive law (recursive least squares)
    %Omegadot = -R*Psi(x(i-1))*(epsilon+e); 
    Rdot     = -R*Psi(x(i-1))*Psi(x(i-1))'*R;
    % % least-squares gradient MRAC adaptive Law
    Omegadot = -Gamma*Psi(x(i-1))*(epsilon + e*P);
    % compute relevant dynamics
    xdot_ref = am*x_ref(i-1) + bm*r;
    xdot     = a *x(i-1) + b*u(i-1) + b*theta*Phi(x(i-1));
    xdot_no_adapt = a*x_no_adapt(i-1) + b*u_no_adapt + b*theta*Phi(x_no_adapt(i-1));
    % integrating/updating 
    x_ref(i) = x_ref(i-1) + xdot_ref*dt;
    x(i)     = x(i-1) + xdot*dt;
    Omega_est(:,i)  = Omega_est(:,i-1) + Omegadot*dt;
    theta_est(i)  = Omega_est(2,i)/b;
    R  = R + Rdot*dt;
    a_est(i)      = Omega_est(1,i);
    x_no_adapt(i) = x_no_adapt(i-1) + xdot_no_adapt*dt;
    e_no_adapt(i) = x_ref(i) - x_no_adapt(i);
    r11(i) = R(1,1);
    r12(i) = R(1,2);
    r21(i) = R(2,1);
    r22(i) = R(2,2);
end

figure('visible','on');
plot(t,x_ref,'b--');
hold on
plot(t,x,'r--');
plot(t,x_no_adapt,'k--');
legend('reference','adapt','no adapt');

figure('visible','on');
plot(t,theta_est);
legend('theta est');

figure('visible','on');
plot(t,a_est,'k');
hold on;
legend('a est');

figure('visible','on');
plot(t,r11,'b');
hold on
plot(t,r12,'k');
plot(t,r21,'r');
plot(t,r22,'m');