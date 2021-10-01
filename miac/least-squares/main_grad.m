% Least squares example 
% Prone to numerical instability  (do not use in real application)
% Joe Chai 2019
Theta = [3; 10; 4; 2]
Phi = @(x) [1; x; x^2; x^3]; % fifth order polynomial
dt = 0.0001; % requires very small time steps or will de-converge (still not robust when using large times)
x = 0:dt:3;
y = zeros(length(x),1);
for i=1:length(x)
    y(i) = Theta'*Phi(x(i));
end

gam = [0.8 0 0 0 ;
       0 0.8 0 0 ;
       0 0 0.8 0 ;
       0 0 0 0.8 ];
Theta_est = zeros(length(Theta),length(x)); Theta_est(:,1) = [1;7;5;0];
e = zeros(length(x),1);
% Compute Params using gradient method
for i=2:length(x)
    e(i-1) = Theta_est(:,i-1)'*Phi(x(i-1)) - y(i-1);
    Thetadot = -gam*Phi(x(i-1))*e(i-1)'  ;
    Theta_est(:,i) = Theta_est(:,i-1) + Thetadot*dt;   
end

figure('visible','on')
plot(x,Theta_est(1,:),'b')
hold on
plot(x,Theta_est(2,:),'k--')
plot(x,Theta_est(3,:),'r')
plot(x,Theta_est(4,:),'m--')
legend('1','2','3','4');

figure('visible','on')
plot(x,e)

Theta_est(:,end)
% This method causes the error to converge but does not really identify the
% params correctly, kind of like mrac