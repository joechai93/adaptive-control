% Demo of direct MRAC for second order siso system
% Dr Nhan Nguyen Text on MRAC section 5.4
% Joe Chai 2019
% Chang line 49 to change between case 1 (A,B unknown) and 2 (A,B known)

clearvars
clear all
% define system
w = 1; zeta = -0.5; b= 1;
A = [0 1;-w*w -2*zeta*w]; B = [0;b]; 
theta = [0.5; -0.1];
phi = @(x) [1;x(1)*x(1)];
Am = [0 1;-2*2 -2*0.5*2]; bm = 4; Bm = [0;bm]; 
% solve lyapunov equation
Q = eye(2);
P = lyap(Am',Q);
Pbar = [P(1,2); P(2,2)];
% sim params
dt = 0.01;
t = 0:dt:50;
x_ref = zeros(2,length(t));
x     = zeros(2,length(t));
u     = zeros(length(t),1);
kx    = zeros(length(t),2);
kr    = zeros(length(t),1);
theta_est  = zeros(2,length(t)); 
e_no_adapt = zeros(2,length(t));
x_no_adapt = zeros(2,length(t));

r = 0;

% adaptive parameters
gam_x = [100 0;0 100];
gam_r = 100;
gam_f = [100 0;0 100];

% non adaptive comparison
kx_true = pinv(B'*B)*B'*(Am - A); %using pseudo inverse for B
kr_true = bm/b;


% run smulation
for i=2:length(t)
    if t(i) > 0.5
        r = sin(2*t(i));
    end
    % control law
    % kx(i-1,:)*x(:,i-1) + kr(i-1)*r - theta_est(:,i-1)'*phi(x(:,i-1)); 
    u(i-1) = kx_true*x(:,i-1) + kr_true*r - theta_est(:,i-1)'*phi(x(:,i-1));
    u_no_adapt = kx_true*x_no_adapt(:,i-1) + kr_true * r;
    % adaptive law
    e = x_ref(:,i-1) - x(:,i-1);
    kxdot = gam_x * x(:,i-1) * e' * Pbar * sign(b);
    krdot = gam_r * r      * e' * Pbar * sign(b);
    theta_estdot  = -gam_f * phi(x(:,i-1)) * e' * Pbar * sign(b); % main one
    % compute relevant dynamics
    xdot_ref = Am*x_ref(:,i-1) + Bm*r;
    xdot     = A*x(:,i-1) + B*u(i-1) + B*theta'*phi(x(:,i-1));
    xdot_no_adapt = A*x_no_adapt(:,i-1) + B*u_no_adapt + B*theta'*phi(x_no_adapt(:,i-1));
    % integrating/updating 
    x_ref(:,i) = x_ref(:,i-1) + xdot_ref*dt;
    x(:,i)     = x(:,i-1) + xdot*dt;
    kx(i,:)    = kx(i-1,:) + (kxdot)'*dt;
    kr(i)    = kr(i-1) + krdot*dt;
    theta_est(:,i)  = theta_est(:,i-1) + theta_estdot*dt;
    x_no_adapt(:,i) = x_no_adapt(:,i-1) + xdot_no_adapt*dt;
    e_no_adapt(:,i) = x_ref(:,i) - x_no_adapt(:,i);
end

figure('visible','on');
plot(t,x_ref(1,:),'b--');
hold on
plot(t,x(1,:),'r--');
plot(t,x_no_adapt(1,:),'k--');
legend('reference','adapt','no adapt');

figure('visible','on');
plot(t,theta_est(1,:));
hold on
plot(t,theta_est(2,:));
legend('theta est');
