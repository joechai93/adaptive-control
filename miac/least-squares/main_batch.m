% Least squares example 
% Joe Chai 2019
Theta = [2; 3; 4; 4; 1; 3]
Phi = @(x) [1; x; x^2; x^3; x^4; x^5]; % fifth order polynomial
x = 0:0.1:5;
y = zeros(length(x),1);
for i=1:length(x)
    y(i) = Theta'*Phi(x(i));
end

% Compute Params using Least Squares
Theta_est = ComputeParamsLS(x,y,Phi)