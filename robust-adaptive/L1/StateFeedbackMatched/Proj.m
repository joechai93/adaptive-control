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

