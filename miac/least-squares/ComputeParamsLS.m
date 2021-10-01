function [Theta,A,B] = ComputeParamsLS(x,y,Phi)
% Chap 6 Dr Nhan Nguyen MRAC
A = 0;
B = 0;
N = length(y);
for i=1:N
    A = A + Phi(x(i))*Phi(x(i))';
    B = B + Phi(x(i))*y(i)';
end
Theta = inv(A)*B;
end
