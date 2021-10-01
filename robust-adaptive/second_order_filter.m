function [yn] = second_order_filter(a0,a1,a2,b2,dt,xn,xn1,xn2,yn1,yn2)
K = 2/dt;
A0 = (b2)/(a0*K*K + a1*K + a2);
A1 = (2*b2)/(a0*K*K + a1*K + a2);
A2 = (b2)/(a0*K*K + a1*K + a2);
B1 = (2*a2 - 2*a0*K*K)/(a0*K*K + a1*K + a2);
B2 = (a0*K*K - a1*K + a2)/(a0*K*K + a1*K + a2);
yn = A0*xn + A1*xn1 + A2*xn2 - B1*yn1 - B2*yn2;
end