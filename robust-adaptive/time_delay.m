% Demo of time delay in MRAC for first order siso system
% Dr Nhan Nguyen Text on MRAC section 8.3
% Includes fix of projection method section 9.2
%
% Joe Chai 2019
%
% This script is supposed to show the effects of time-delay on robustness
% of MRAC. Instability can result when adaptive parameter exceeds a limit.

clearvars
clear all
% define system
a = -1; b = 1; td = 0.5; 
am = -2; bm = 2;

% sim params
dt = 0.01;
t = 0:dt:10;
t_buff = td/dt;
x     = zeros(length(t),1); x(1) = 1;
u     = zeros(length(t),1);
kx    = zeros(length(t),1); kx(1) = (am - a)/b;
w     = zeros(length(t),1);
r = 0;

% adaptive parameters
gam_x = 2.9; % 2 still stable, 2.9 unstable

% projection method Example 9.2
g = @(k) abs(k - 0) - 3;
use_projection_method = true; % toggle between true and false to test

% run smulation
for i=2:length(t)
    % control law
    u(i-1) = kx(i-1)*x(i-1);
    % adaptive law
    if (g(kx(i-1)) < 0 || (abs(g(kx(i-1)) - 0)<0.001 && -x(i-1)*x(i-1)*b*1 <= 0)) ...
        || use_projection_method==false
    
        kxdot = -gam_x * x(i-1) * x(i-1) * b;
    else
        kxdot = 0; % del(g_k) = 1 so the projection is zero
    end
    % compute relevant dynamics
    if t(i) > td
        xdot     = a*x(i-1) + b*u(i-1-t_buff);
    else
        xdot = 0;
    end
    % integrating/updating 
    x(i)     = x(i-1) + xdot*dt;
    kx(i)    = kx(i-1) + kxdot*dt;
end

figure('visible','on');
subplot(2,2,1)
plot(t,x);
legend('x')
subplot(2,2,2)
plot(t,kx);
legend('kx')
subplot(2,2,4)
plot(t,u)
legend('u')
