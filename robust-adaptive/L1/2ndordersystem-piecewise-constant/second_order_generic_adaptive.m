clearvars

A = [Mq Ma;
     1 -Na/V];
Bm = [Mde;
      -Nde/V];
  
Bum = [0.0005; -0.0005*Bm(1)/Bm(2)];
B = [Bm Bum];


% efficient model predictive control
[Vec,Val] = eig(A);
lambda1 = Val(1,1);
lambda2 = Val(2,2);
S = Vec;
dt_mpc = 0.5;
STM = real(S*[exp(lambda1*dt_mpc) 0;
        0   exp(lambda2*dt_mpc)]*inv(S));
C = [0 1];
Q = 0.99;%0.99 0;0 0.99];
R = 0.001;
F = C*STM;
I = [1 0; 0 1];
G = C*inv(A)*(STM-I)*Bm;
K = inv(G'*Q*G + R) * G' * Q;

Am = A - Bm*K*F;

n = 2;
m = 1;

Ad = zeros(2,2);%[-0.8*Mq 0.8*Ma;
     %  0.3 1.0];

w0 = 1.0;%0.4;


%Kg = -inv(C*inv(Am)*Bm);

% sim params
dt = 0.001;
t = 0:dt:40;
t_delay = 0.0; t_delay_tick = 0;
x_no_adapt  = zeros(n,length(t)); x_no_adapt(2,2) = 2; x_no_adapt(2,1) = 2;

u_no_adapt  = zeros(m,length(t));
u_no_adapt_delay = zeros(m,length(t));

r           = zeros(m,length(t));

% reference ramp defines
tstart = 15;
tmid = 20;
tmid2 = 30;
tend = 35;
rhigh = 5;
rlow = 0;

for i=3:length(t)
     %r(i) = sin(2*t(i)/pi);
    % reference ramp
    if t(i) < tstart
        r(i) = 0;
    elseif t(i) <= tmid && t(i) > tstart
        r(i) = ((t(i)-tstart)/(tmid-tstart))*(rhigh-rlow) + rlow;
    elseif t(i) > tmid && t(i) <= tmid2
        r(i) = rhigh;
    elseif t(i) > tmid2 && t(i) < tend
        r(i) = ((t(i)-tmid2)/(tend-tmid2))*(rlow-rhigh) + rhigh;
    elseif t(i) > tend
        r(i) = rlow;
    end
    
    u_no_adapt(:,i-1) = K*(r(i)-F*x_no_adapt(:,i-1));

    % tdm analysis
    if t_delay_tick > t_delay/dt
        u_no_adapt_delay(i-1) = u_no_adapt(i-t_delay_tick);
    else
        t_delay_tick = t_delay_tick + 1;
    end
    
    % compute relevant dynamics
   
    % switch for tdm analysis
    xdot_no_adapt  = (A + Ad)*x_no_adapt(:,i-1) + Bm*w0*u_no_adapt_delay(:,i-1);      
    
    % integrating/updating 
    x_no_adapt(:,i)  = x_no_adapt(:,i-1) + xdot_no_adapt*dt;
     
end

figure('visible','on');
plot(t, r, 'k-');
hold on
plot(t, x_no_adapt(2,:), 'g--');
grid on;
legend('r','y');

figure('visible','on');
plot(t,u_no_adapt_delay);
grid on;
legend('u delay');