%L1NORM
%
%   function [L1norm,err,U,L,n,tol]=l1norm(G,tol);
%
%	l1norm(G)
%
% 	Calculate L1-norm (Rutland & Lane's algorithm)
% 	of the impulse response of a proper SISO LTI system with
%	minimal state space realization G=(A,B,C,0)
%
%	[norm,err,U,L,n,tol]=l1norm(G)
%
%	Additional outputs:
%	err - error estimate
%	U - upper bound on norm
%	L - lower bound on norm
%	niter - number of iterations
%	tol = -alpha * 1e-3 where alpha is maximum
%	      real part of the eigenvalues of A
%
%	l1norm(G,tol)
%
% 	Calculate L1-norm by Rutland & Lane's algorithm
% 	of the impulse response of a SISO LTI system with
%	minimal state space realization (A,B,C) to a
%	tolerance tol
%   (default tol=-alpha/1000 where alpha is the abscissa of stability)
%
%	l1norm(G,tol,maxiter)
%
% 	Calculate L1-norm by Rutland & Lane's algorithm
% 	of the impulse response of a SISO LTI system with
%	minimal state space realization (A,B,C) to a
%	tolerance tol with a maximum number of of iterations niter 
%   (default maxiter=20) 
%   NOTE that memory problems can occur for niter>20, so increase with
%        caution
%   
%
%	See Rutland N.K. & Lane P.G, "Computing the 1-norm of the impulse 
%        response of linear time-invariant systems", 
%        Systems and Control Letters, Volume 26, Number 3, pp. 211-221 
%        http://dx.doi.org/10.1016/0167-6911(95)00022-2
% 

% (c) JF Whidborne 28/4/95 
%     modified jfw 1/5/2013, 6/10/14, 13/11/14
%     

function [L1norm,err,U,L,tol,niter]=l1norm(G,tol,maxiter)
narginchk(1,3);

Gss=ss(G); % ensure state space realization

% check SISO
[ny,nu] = size(Gss.d);
if ny~=1||nu~=1
    error('Error: System not SISO')
end

% Check strict properness
if Gss.d~=0
    error('Error: System not strictly proper')
end

A=Gss.a;B=Gss.b;C=Gss.c; % get state matrices/vectors
% abscissa of stability alpha
alpha=max(real(eig(A)));
if alpha>0, error('unstable system');end
if cond(A)>1e6, % 13/11/14
    warning('ill conditioned system - results may be inaccurate - use balreal to obtain a balanced realization');
end

if nargin<3, maxiter=20; end% maximum iterations
if nargin<2, tol=-alpha*.001; end% tolerance
%OPTIONS1=[0,-alpha*1e-4];
[ns,~] = size(A); % state dimension
x0 = zeros(1,ns);% initial condition
G0=-C*inv(A)*B;% final steady state value (step response)

% 1st step T=0, N=0
%[sig_min,OPTIONS]=fmin('normtail',0,-alpha,OPTIONS1,A,B',C,ns);

[~,nu_u_Nplus1] = fminbnd(@normtail,0,-alpha,...
    optimset('Display','on','TolFun',-alpha*1e-4),A,B',C,ns);
%fplot(@(sig) normtail(sig,A,B',C,ns),[0,-alpha]);

nu_u(1)= nu_u_Nplus1;%  weighted H2 norm
U=nu_u(1); % upper bound
nu_l(1)=abs(G0); % steady state for step response
L=nu_l(1); % lower bound

niter=1;N=1;T=-5/alpha; T_flag=1; % estimate of required time T
%tolerance_error=.5*(U-L)/L;

while((((U-L)/L/2)>=tol) && niter<maxiter)% increase maxiter for accuracy if memory available !!
    niter=niter+1; % iteration number
    h=T/N; % step length
    
    [Phi, Gamma] = c2d(A,B,h); % get difference equation
    
    % evaluate nu_l(n) lower bound
    
    y = ltitr(Phi,Gamma,[1;zeros(N,1)],x0)*C.';% get impulse response
    nu_l_k=abs(y(2:N+1));% lower bounds in intervals
    if T_flag, % time has changed - get nu_l_{N+1} error
        %y1 = ltitr(Phi,Gamma,ones(N+1,1),x0)*C.';% get 1 response
        %nu_l_Nplus1= abs(G0- y1(N+1));
        nu_l_Nplus1= abs(G0- sum(y));%sum(y) is step response at t=T
    end% if
    sum_nu_l_k=sum(nu_l_k); % lower bound sum <=T
    nu_l(niter) = sum_nu_l_k + nu_l_Nplus1;% eqn(8)
    
    % evaluate nu_u(n) upper bound
    
    x = ltitr(Phi,Gamma,zeros(N+1,1),B);% state impulse response, get value of last point
    CTC=C'*C; W=lyap(A',CTC-Phi'*CTC*Phi);% Lyapanov eqn(36)
    %nu_u_k=(h*sum(((x(1:N,:)*W).*x(1:N,:))')).^.5;% unrefined upper bounds
    %sum_nu_u_k=sum(nu_u_k); % unrefined upper bound sum <=T
    % now evaluate refined upper bounds
    nu_u_k=nu_l_k;
    AWA=A'*W*A;
    int_diff=sqrt(h*sum(((x(1:N,:)*AWA).*x(1:N,:))'));% eqn (41)
    e_abs=abs(C*x');
    %max_e=max(e_abs(1:N),e_abs(2:N+1));
    n_rub=find((max(e_abs(1:N),e_abs(2:N+1)))<=int_diff);% intervals which do not have exact lower bound
    nu_u_k(n_rub)=sqrt(h*sum(((x(n_rub,:)*W).*x(n_rub,:))'));% refined upper bounds
    
    if T_flag, % time has changed - get nu_u_{N+1} error
%        [sig_min,OPTIONS]=fmin('normtail',0,-alpha,OPTIONS1,A,x(N+1,:),C,ns);
        [~, nu_u_Nplus1] = fminbnd(@normtail,0,-alpha,...
            optimset('Display','on','TolFun',-alpha*1e-4),A,x(N+1,:),C,ns);
        T_flag=0; % reset flag
    end% if
    sum_nu_u_k=sum(nu_u_k); % upper bound sum <=T
    nu_u(niter) = sum_nu_u_k + nu_u_Nplus1;% eqn(11)
    
    % calculate lower bound L
    L=max(nu_l);
    % calculate upper bound U
    U=min(nu_u);
    
    N=2*N;
    %double T criterion - change final time and set flag
    if((sum_nu_u_k-sum_nu_l_k)<(nu_u_Nplus1-nu_l_Nplus1)),T=2*T;T_flag=1;end
end%while
err=(U-L)/L/2;% error estimate
L1norm=(U+L)/2; % mean of upper & lower bounds
return


function nu=normtail(sig,A,x,C,ns)
% minimising function used by l1norm.m - eqn(39)
%
% 1/2\sigma \int_0^\infty |\exp(\sigma t) e(t+T,\delta)|^2 dt
%

W=lyap(A'+sig*eye(ns),C'*C);
nu= ((x*W*x')/2/sig)^(0.5);%need square root
return
