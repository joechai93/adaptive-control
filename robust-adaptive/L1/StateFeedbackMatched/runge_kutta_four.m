function [xout] = runge_kutta_four(xin, xdot, u, dt)
% runge kutta integrator for state space model
% xin: current state vector
% xdot: f(x,u) state dynamics function
% u: input vector
% dt: integration time step

xdot1 = xdot(xin,u);
An = xdot1*dt;

xdot2 = xdot(xin + 0.5*An, u);
Bn = xdot2*dt;

xdot3 = xdot(xin + 0.5*Bn, u);
Cn = xdot3*dt;

xdot4 = xdot(xin + Cn, u);
Dn = xdot4*dt;

xout = xin + (An + 2.0*Bn + 2.0*Cn + Dn)/6.0;
end

%%void VehicleErrorModel::Integrate(Eigen::VectorXd& X, Eigen::VectorXd U, double DT)
%%{
%%    Eigen::VectorXd Xdot1 = StateRates(X, U);	
%%    Eigen::VectorXd An = Xdot1 * DT;
%%
%%    Eigen::VectorXd Xdot2 = StateRates(X + 0.5*An, U);
%%    Eigen::VectorXd Bn = Xdot2 * DT;
%%
%%    Eigen::VectorXd Xdot3 = StateRates(X + 0.5*Bn, U);
%%    Eigen::VectorXd Cn = Xdot3 * DT;
%%
%%    Eigen::VectorXd Xdot4 = StateRates(X + Cn, U);
%%    Eigen::VectorXd Dn = Xdot4 * DT;
%%
%%    X = X + (An + 2.0*Bn + 2.0*Cn + Dn)/6.0;
%%
%%}

