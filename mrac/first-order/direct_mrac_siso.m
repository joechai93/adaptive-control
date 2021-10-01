% Demo of direct MRAC for first order siso system
% Dr Nhan Nguyen Text on MRAC section 5.2
% Joe Chai 2019
% Assumed a,b known but theta unknown (Case 2)
clearvars
clear all
% define system
a = 4; b = 2; theta = 0.8;
am = -2; bm = 2;

% sim params
dt = 0.01;
t = 0:dt:50;
x_ref = zeros(length(t),1);
x     = zeros(length(t),1);
u     = zeros(length(t),1);
kx    = zeros(length(t),1);
kr    = zeros(length(t),1);
theta_est  = zeros(length(t),1);
e_no_adapt = zeros(length(t),1);
x_no_adapt = zeros(length(t),1);

r = 0;

% adaptive parameters
gam_x = 0.1;
gam_r = 0.1;
gam_f = 0.6;

% non adaptive comparison
kx_true = (am - a)/b;
kr_true = bm/b;

% run smulation
for i=2:length(t)
    if t(i) > 0.5
        r = sin(t(i));
    end
    % control law
    u(i-1) = kx_true*x(i-1) + kr_true*r - theta_est(i-1)*x(i-1); %kx(i-1)*x(i-1) + kr(i-1)*r - fest(i-1)*x(i-1); method one
    u_no_adapt = kx_true*x_no_adapt(i-1) + kr_true * r;
    % adaptive law
    e = x_ref(i-1) - x(i-1);
    kxdot = gam_x * x(i-1) * e * sign(b);
    krdot = gam_r * r      * e * sign(b);
    theta_estdot  = -gam_f * x(i-1) * e * sign(b); % main one
    % compute relevant dynamics
    xdot_ref = am*x_ref(i-1) + bm*r;
    xdot     = a *x(i-1) + b*u(i-1) + b*theta*x(i-1);
    xdot_no_adapt = a*x_no_adapt(i-1) + b*u_no_adapt + b*theta*x_no_adapt(i-1);
    % integrating/updating 
    x_ref(i) = x_ref(i-1) + xdot_ref*dt;
    x(i)     = x(i-1) + xdot*dt;
    kx(i)    = kx(i-1) + kxdot*dt;
    kr(i)    = kr(i-1) + krdot*dt;
    theta_est(i)  = theta_est(i-1) + theta_estdot*dt;
    x_no_adapt(i) = x_no_adapt(i-1) + xdot_no_adapt*dt;
    e_no_adapt(i) = x_ref(i) - x_no_adapt(i);
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
