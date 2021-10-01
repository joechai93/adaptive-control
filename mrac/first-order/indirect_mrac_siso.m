% Demo of indirect MRAC for first order siso system
% Dr Nhan Nguyen Text on MRAC section 5.3
% Joe Chai 2019
clearvars
close all
% define system
a = 4; b = 2; theta = 0.8;
am = -2; bm = 2;

% sim params
dt = 0.01;
t = 0:dt:50;
x_ref = zeros(length(t),1);
x     = zeros(length(t),1);
u     = zeros(length(t),1);
a_est      = zeros(length(t),1); a_est(1) = a;
b_est      = zeros(length(t),1); b_est(1) = b;
theta_est  = zeros(length(t),1); theta_est(1) = 0.1;
e_no_adapt = zeros(length(t),1);
x_no_adapt = zeros(length(t),1);

r = 0;

% adaptive parameters
gam_a = 0.1;
gam_b = 0.1;
gam_f = 0.6;
b0 = 0.1;

% true gains
kx_true = (am - a)/b;
kr_true = bm/b;

% run smulation
for i=2:length(t)
    if t(i) > 0.5
        r = sin(t(i));
    end
    % control law
    kx = (am - a_est(i-1))/b_est(i-1);
    kr = bm/b_est(i-1);
    u(i-1) = kx*x(i-1) + kr*r - theta_est(i-1)*x(i-1); 
    u_no_adapt = kx_true*x_no_adapt(i-1) + kr_true * r;
    % adaptive law
    e = x_ref(i-1) - x(i-1);
    adot = -gam_a * x(i-1) * e;
    if (b_est(i-1) >= b0) % Eqn 5.54
        bdot = -gam_b * u(i-1) * e;
        if (b_est(i-1) == b0 && bdot < 0)
            bdot = 0;
        end
    end
    theta_estdot  = -gam_f * x(i-1) * e * sign(b); % main one
    % compute relevant dynamics
    xdot_ref = am*x_ref(i-1) + bm*r;
    xdot     = a *x(i-1) + b*u(i-1) + b*theta*x(i-1);
    xdot_no_adapt = a*x_no_adapt(i-1) + b*u_no_adapt + b*theta*x_no_adapt(i-1);
    % integrating/updating 
    x_ref(i) = x_ref(i-1) + xdot_ref*dt;
    x(i)     = x(i-1) + xdot*dt;
    a_est(i)    = a_est(i-1) + adot*dt;
    b_est(i)    = b_est(i-1) + bdot*dt;
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

figure('visible','on');
plot(t,a_est,'k');
hold on;
plot(t,b_est,'r');
legend('a est','b est');