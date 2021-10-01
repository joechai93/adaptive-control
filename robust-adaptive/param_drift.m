% Demo of parameter drift in MRAC for first order siso system
% Dr Nhan Nguyen Text on MRAC section 8.1
% Includes fix using dead zone method section 9.1
%
% Joe Chai 2019
%
% This script is supposed to show the parameter drift
% phenomenon of adaptive control under enough unmodelled disturbance
% Kx is unbounded and can drift indefinitely, leading to eventual instability
clearvars
clear all
% define system
% regulator system (ref model is zeros)
a = 4; b = 2; 
am = -2; bm = 2;

% sim params
dt = 0.01;
t = 0:dt:5000;
x     = zeros(length(t),1); x(1) = 1;
u     = zeros(length(t),1);
kx    = zeros(length(t),1); kx(1) = (am - a)/b;
w     = zeros(length(t),1);
r = 0;

% adaptive parameters
gam_x = 10;

% dead zone method
% sec 9,1
dead_zone_on = false;
e0 = 0.1;

% run smulation
for i=2:length(t)
    w(i-1) = 3.4*exp(-0.0001*(t(i)));
    % control law
    u(i-1) = kx(i-1)*x(i-1);
    % adaptive law
    if abs(x(i-1)) > e0 || dead_zone_on == false
        kxdot = -gam_x * x(i-1) * x(i-1) * b;
    else
        kxdot = 0;
    end
    % compute relevant dynamics
    xdot     = a *x(i-1) + b*u(i-1) + w(i-1);
    % integrating/updating 
    x(i)     = x(i-1) + xdot*dt;
    kx(i)    = kx(i-1) + kxdot*dt;
end

figure('visible','on');
subplot(2,2,1)
plot(t,x,'k--');
legend('x')
subplot(2,2,2)
plot(t,kx);
legend('kx')
subplot(2,2,3)
plot(t,w)
legend('w')
subplot(2,2,4)
plot(t,u)
legend('u')
