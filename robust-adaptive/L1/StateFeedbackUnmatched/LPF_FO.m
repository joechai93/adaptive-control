function [yn] = LPF_FO(un,y,dt,wc)
yn = (dt*wc/(1+dt*wc))*un + (1/(1+dt*wc))*y;
end
