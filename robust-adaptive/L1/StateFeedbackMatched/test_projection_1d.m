% Testing aspects of 1d (non zero centred) projection method 
% in adaptive control
% mainly for uncertain input gain
% Joe Chai May 2019
clearvars

dOm = 0.05;
Omega_ = 0.4:dOm:1.9;
Omega0_ = [0.5 1.8];

Omega_mid = (Omega0_(1) + Omega0_(2))/2;
Omega_d_max = abs(Omega0_(2) - Omega_mid);

f_ = zeros(length(Omega_),1);
delf_ = zeros(length(Omega_),1);
proj_ = zeros(length(Omega_),1);


for i=1:length(Omega_)
    omega = Omega_(i)-Omega_mid;
    f_(i) = f(omega,Omega_d_max);
    delf_(i) = delf(omega,Omega_d_max);
    proj_(i) = Proj(omega,1.0,Omega_d_max); 
end

figure()
plot(Omega_,proj_);
legend('unit projection');
%quiver(theta1_,-theta1_,proj1_',proj2_');

figure()
plot(Omega_,f_);
hold on
plot(Omega_,delf_);
legend('f','delf');

function [fout] = f(Theta,ThetaMax)
% convex compact set 
% ThetaMax must be scalar - norm bound imposed on Theta
assert(length(ThetaMax)==1);
e0 = 0.1;
fout = ((e0 + 1)*Theta'*Theta - ThetaMax*ThetaMax)/(e0*ThetaMax*ThetaMax);
end

function [delfout] = delf(Theta,ThetaMax)
e0 = 0.1;
assert(length(ThetaMax)==1);
delfout = ((2*(e0+1))/(e0*ThetaMax*ThetaMax))*Theta;
end

function [proj_out] = Proj(Theta,Y,ThetaMax)
delf_now = delf(Theta,ThetaMax);
% pg 292 L1 Adaptive Control Text by Dr Hovakimyan
if f(Theta,ThetaMax) < 0
    % Inside convex set, no modification 
    ret = Y;   
elseif f(Theta,ThetaMax) >= 0 && delf_now'*Y <= 0
    % outside convex set but gradient is bringing it back
    ret = Y;
else
    % outside convex set and gradiant is pushing it further
    ret = Y - ((delf_now)*delf_now'/(norm(delf_now)^2))*Y*f(Theta,ThetaMax);
end
proj_out = ret;
end