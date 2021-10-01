% Simulation of rohrs example
% demonstrating MRAC lack of robustness to unmodelled dynamics
% 2.3.4 Dr Hovakimyan Text
addpath('./L1/StateFeedbackMatched');

r1 = @(t) 0.3 + 1.85*sin(16.1*t);
r2 = @(t) 0.3 + 2*sin(8*t);
r = r1;

dt = 0.00001;
t = 0:dt:20;
y = zeros(length(t),1);
ym = zeros(length(t),1);
mu = zeros(length(t),1);
u = zeros(length(t),1);
ky = zeros(length(t),1); ky(1) = -0.65; ky(2) = -0.65; 
kr = zeros(length(t),1); kr(1) = 1.14; kr(2) = 1.14; 


for i=3:length(t)
    % control law
    u(i)  = ky(i-1)*y(i-1) + kr(i-1)*r(t(i-1));
    % adaptive law
    e = y(i-1) - ym(i-1);
    kydot = -e*y(i-1);
    krdot = -e*r(t(i-1));
    % time marching
    mu(i) = second_order_filter(1,30,229,229,dt,u(i),u(i-1),u(i-2),mu(i-1),mu(i-2));
    y(i)  = second_order_filter(0,1,1,2,dt,mu(i),mu(i-1),mu(i-2),y(i-1),y(i-2));
    %y(i)  = second_order_filter(0,1,1,2,dt,u(i),u(i-1),u(i-2),y(i-1),y(i-2));
    ym(i) = LPF_FO(r(t(i)),ym(i-1),dt,3);
    ky(i) = ky(i-1) + kydot*dt;
    kr(i) = kr(i-1) + krdot*dt;
end

figure()
subplot(2,1,1);
plot(t,ym);
hold on
plot(t,y);
grid on;
legend('ym','y');

subplot(2,1,2);
plot(t,ky);
hold on
plot(t,kr);
grid on;
legend('ky','kr');

