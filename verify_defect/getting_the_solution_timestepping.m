function [a, d_N, d_infty, a_end] = getting_the_solution_timestepping(N,n,tspan,v_hat,rigorous)

% theta = angle/180*pi;
% gamma = exp(1i*mid(theta));
gamma = 1i;
stepsize = tspan(2)-tspan(1);
%
% N: # of fourier mode
% tmax: end of time

chebfunpref.setDefaults('factory');
chebfunpref.setDefaults('fixedLength',n);
%

% opts = odeset('abstol',1e-18,'reltol',2.22045e-14);
opts = odeset('abstol',1e-18,'reltol',1e-18);

u_cheb = chebfun.ode45(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
a = chebcoeffs(u_cheb);
a_end = u_cheb(end);
% u_cheb = chebfun.ode113(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
% u_cheb = chebfun.ode15s(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);


% k = (-N:N);
% weight=1.05;
% nu = (weight).^(-abs(k));
% v_hat = v_hat./nu;
% u_cheb = chebfun.ode45(@(t,y) ode_func2(y,gamma,nu'),tspan,v_hat,opts);
% a = chebcoeffs(u_cheb).*nu;
% a_end = u_cheb(end).*nu;

chebfunpref.setDefaults('factory');

rescaleFactork = stepsize/2;
du = ChebDerCoeffs(a,rigorous)/rescaleFactork;
du = [du;zeros(1,size(du,2))];

du = [du(1,:);du(2:end,:)/2];

a0 = [flipud(a(2:end,:))/2;a(1,:);a(2:end,:)/2];
% a0(:,1:N) = fliplr(a0(:,N+2:end));
a0(:,N+2:end) = fliplr(a0(:,1:N)); % symmetry-condition


if rigorous>0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evaluattion of the residual with interval arithmetic %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Convolution %%%
  a = a0; ia = intval(a); ipi = intval('pi');
%   itheta = theta;
%   itheta = ipi*angle/180; 
%   igamma = exp(1i*itheta);
  igamma = 1i;
  [~,ia2full] = convtensor(ia,ia);
  
  ia = ia(n:end,:);
  ia2full = ia2full(2*n-1:end,:);
  
  ia_full = [zeros(n,N) ia zeros(n,N);zeros(n-1,4*N+1)];
  du_full = [zeros(n,N) du zeros(n,N);zeros(n-1,4*N+1)];
  k_full = -2*N:2*N;
  ires_full = du_full + igamma*((4*ipi^2)*(k_full.^2).*ia_full-ia2full);

%   nu_full =  (weight).^(abs(k_full)); % weighted norm test
%   ires_full = ires_full.*nu_full;
  
  ires_full(2:end,:) = 2*ires_full(2:end,:);% back to two-sided Chebyshev coefficients
  
  %   disp(sup(norm(ires_full,1)));
  
  d_N = mag(sum(sum(abs(ires_full(:,N+1:3*N+1)))));
  ires_full(:,N+1:3*N+1)=0;
  d_infty = mag(sum(sum(abs(ires_full))));
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Evaluation of the residual without interval arithmetic %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Convolution %%%
  a = a0;  
  [~,a2full] = convtensor(a,a);
  
  a = a(n:end,:);
%   a2 = a2(n:end,:);
  a2full = a2full(2*n-1:end,:);
  
%   k = -N:N;
%   res = du + gamma*((4*pi^2)*(k.^2).*a-a2);
  a_full = [zeros(n,N) a zeros(n,N);zeros(n-1,4*N+1)];
  du_full = [zeros(n,N) du zeros(n,N);zeros(n-1,4*N+1)];
  k_full = -2*N:2*N;
  res_full = du_full + gamma*((4*pi^2)*(k_full.^2).*a_full-a2full);
  
%   res(2:end,:) = 2*res(2:end,:);% back to chebyshev coefficients
  res_full(2:end,:) = 2*res_full(2:end,:);% back to chebyshev coefficients
%   disp(sum(sum(abs(res),2)))
%   disp(sum(sum(abs(res_full),2)))
%   disp(sum(sum(abs(res_full),2))-sum(sum(abs(res),2)))
  d_N = sum(sum(abs(res_full(:,N+1:3*N+1))));
  res_full(:,N+1:3*N+1)=0;
  d_infty = sum(sum(abs(res_full)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% We output the data one-sided in Chebyshev %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = a0(n:end,:);

chebfunpref.setDefaults('factory');
end

function dy = ode_func(y,gamma)

N = size(y,1);
fy = quadratic(y,y);
m = (N-1)/2;
k = (-m:m)';
dy = -(4*pi^2)*(k.^2).*y+fy;
dy = dy*gamma;

end
% function dy = ode_func2(y,gamma,nu)
% 
% N = size(y,1);
% m = (N-1)/2;
% k = (-m:m)';
% % nu = (omega).^(-abs(k));
% fy = quadratic(nu.*y,nu.*y);
% % k = (-m:m)';
% dy = -(4*pi^2)*(k.^2).*y+fy./nu;
% dy = dy*gamma;
% 
% end