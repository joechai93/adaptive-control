% Demo of indirect MRAC for second order siso system
% Dr Nhan Nguyen Text on MRAC section 5.5
% Joe Chai 2019

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
w_est      = zeros(length(t),1); w_est(1) = 1;
zeta_est   = zeros(length(t),1); zeta_est(1) = -0.5;
b_est      = zeros(length(t),1); b_est(1) = 1;
theta_est  = zeros(2,length(t)); %theta_est(:,1) = theta;
e_no_adapt = zeros(2,length(t));
x_no_adapt = zeros(2,length(t));

r = 0;

% adaptive parameters
gam_w = 1;
gam_z = 1;
gam_b = 1;
gam_f = [30 0;0 30];
w0 = 0.02;
b0 = 0.1;
% non adaptive comparison
kx_true = pinv(B'*B)*B'*(Am - A); %using pseudo inverse for B
kr_true = bm/b;


% run simulation
for i=2:length(t)
    if t(i) > 0.5
        r = sin(2*t(i));
    end
    % control law
    B_est = [0; b_est(i-1)];
    A_est = [0 1;-w_est(i-1)*w_est(i-1) -2*zeta_est(i-1)*w_est(i-1)];
    kx = pinv(B_est'*B_est)*B_est'*(Am - A_est);
    kr = bm/b_est(i-1);
    % kx(i-1,:)*x(:,i-1) + kr(i-1)*r - theta_est(:,i-1)'*phi(x(:,i-1)); 
    u(i-1) = kx*x(:,i-1) + kr*r - theta_est(:,i-1)'*phi(x(:,i-1));
    u_no_adapt = kx_true*x_no_adapt(:,i-1) + kr_true * r;
    % adaptive law
    ubar = kx*x(:,i-1) + kr*r;
    e = x_ref(:,i-1) - x(:,i-1);
    if w_est(i-1) >= w0 % Eqn 5.126
        wdot = gam_w * x(1,i-1) * e' * Pbar / (2*w_est(i-1));
        if (w_est(i-1)==w0 && wdot < 0)
            wdot = 0;
        end
    else
        wdot = 0;
    end
    zetadot = (gam_z*x(2,i-1)*w_est(i-1) - gam_w*x(1,i-1)*zeta_est(i-1)) * e' * Pbar / (2*w_est(i-1)*w_est(i-1));
    if (b_est(i-1) >= b0) 
        bdot  = -gam_b * ubar * e' * Pbar;
        if (b_est(i-1) == b0 && bdot < 0)
            bdot = 0;
        end
    else
        bdot = 0;
    end
    theta_estdot  = -gam_f * phi(x(:,i-1)) * e' * Pbar * sign(b); % main one
    % compute relevant dynamics
    xdot_ref = Am*x_ref(:,i-1) + Bm*r;
    xdot     = A*x(:,i-1) + B*u(i-1) + B*theta'*phi(x(:,i-1));
    xdot_no_adapt = A*x_no_adapt(:,i-1) + B*u_no_adapt + B*theta'*phi(x_no_adapt(:,i-1));
    % integrating/updating 
    x_ref(:,i) = x_ref(:,i-1) + xdot_ref*dt;
    x(:,i)     = x(:,i-1) + xdot*dt;
    w_est(i,:)    = w_est(i-1,:) + wdot*dt;
    zeta_est(i)    = zeta_est(i-1) + zetadot*dt;
    b_est(i)       = b_est(i-1) + bdot*dt;
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

figure('visible','on');
plot(t,w_est,'b')
hold on
plot(t,b_est,'k')
plot(t,zeta_est,'r')
legend('w','b','zeta');