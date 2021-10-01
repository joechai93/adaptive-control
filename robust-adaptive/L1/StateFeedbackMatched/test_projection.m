% Testing aspects of projection method 
% in adaptive control
% Joe Chai May 2019
clearvars



dTh = 0.05;
theta1_ = -0.9:dTh:0.9;
theta2_ = -0.9:dTh:0.9;
theta_max = 1;
f_ = zeros(length(theta1_),length(theta2_));
delf1_ = zeros(length(theta1_)*length(theta2_),1);
delf2_ = zeros(length(theta1_)*length(theta2_),1);
thetaX_ = zeros(length(theta1_)*length(theta2_),1);
thetaY_ = zeros(length(theta1_)*length(theta2_),1);
fplt_ = zeros(length(theta1_)*length(theta2_),1);
proj1_ = zeros(length(theta1_),1);
proj2_ = zeros(length(theta1_),1);

for i=1:length(theta1_)
    for j=1:length(theta2_)
        theta = [theta1_(i); theta2_(j)];
        f_(i,j) = f(theta,theta_max);
        thetaX_((i-1)*length(theta1_) + j) = theta1_(i);
        thetaY_((i-1)*length(theta2_) + j) = theta2_(j);
        fplt_((i-1)*length(theta2_) + j) = f_(i,j);
        delfnow = delf(theta,theta_max);
        delf1_((i-1)*length(theta2_) + j) = delfnow(1);
        delf2_((i-1)*length(theta2_) + j) = delfnow(2);
    end
    proj_now = Proj([theta1_(i);-theta1_(i)],[0.2;0.5],theta_max); 
    proj1_(i) = proj_now(1);
    proj2_(i) = proj_now(2);
end

figure()
plot3(thetaX_,thetaY_,fplt_,'o')

figure()
contour(theta1_,theta2_,f_);
hold on
%quiver(thetaX_,thetaY_,delf1_,delf2_);
quiver(theta1_,-theta1_,proj1_',proj2_');

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
    ret = Y - ((delf_now)*delf_now'/(norm(delf_now)^2))*Y;%*f(Theta,ThetaMax);
end
proj_out = ret;
end