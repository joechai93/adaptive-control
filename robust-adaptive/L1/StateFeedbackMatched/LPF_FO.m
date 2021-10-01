function [yn] = LPF_FO(un,y,dt,wc)
wc  = wc * 6.28; % hz to rad/s
yn = (dt*wc/(1+dt*wc))*un + (1/(1+dt*wc))*y;
end
