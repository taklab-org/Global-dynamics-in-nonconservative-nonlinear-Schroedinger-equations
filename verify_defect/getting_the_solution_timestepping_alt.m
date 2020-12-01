function [a, d_0, d_tail, a_end] = getting_the_solution_timestepping_alt(N,n,tspan,v_hat,rigorous,res_iter)
arguments
  N;  n;  tspan;  v_hat;
  rigorous = true;
  res_iter = false;
end

gamma = 1i;
h = tspan(2)-tspan(1);
%
% N: # of fourier mode
% tmax: end of time

% Compute approximation by chebfun
chebfunpref.setDefaults('factory');
chebfunpref.setDefaults('fixedLength',n);
% opts = odeset('abstol',1e-18,'reltol',2.22045e-14);
opts = odeset('abstol',1e-18,'reltol',1e-18);
u_cheb = chebfun.ode45(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
% u_cheb = chebfun.ode113(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
% u_cheb = chebfun.ode15s(@(t,y) ode_func(y,gamma),tspan,v_hat,opts);
% options = bvpset('abstol',1e-18,'RelTol',1e-18);
% u_cheb = bvp4c(@(t,y) ode_func(y,gamma),@(va,vb) va(:)-v_hat.',u_cheb,options);
% u_cheb = bvp5c(@(t,y) ode_func(y,gamma),@(va,vb) va(:)-v_hat.',u_cheb,options);

chebfunpref.setDefaults('factory');

%%%%

[a0, d_0, d_tail, a_end] = residual_estimation(u_cheb,N,n,h,rigorous);
%   disp(['delta_0 = ',num2str(d_0)])
%   disp(['delta_tail = ',num2str(d_tail)])
% 
if res_iter
  for iter_num = 1:3
    chebfunpref.setDefaults('fixedLength',n);
%     opts = odeset('abstol',1e-17,'reltol',1e-18);
%     u_collec = chebfun.ode113(@(t,y) ode_res_iter(t,y,u_cheb,gamma),tspan,zeros(1,2*N+1),opts);
    u_collec = chebfun.ode15s(@(t,y) ode_res_iter(t,y,u_cheb,gamma),tspan,zeros(1,2*N+1),opts);
%     u_collec = chebfun.ode45(@(t,y) ode_res_iter(t,y,u_cheb,gamma),tspan,zeros(1,2*N+1),opts);
%     options = bvpset('abstol',1e-16,'RelTol',1e-18);
%     options = bvpset('RelTol',1e-18);
%     u_collec = bvp4c(@(t,y) ode_res_iter(t,y,u_cheb,gamma),@(va,vb) va(:),0*u_cheb,options);
    u_cheb = u_cheb - u_collec;
    chebfunpref.setDefaults('factory');
    [a0, d_0, d_tail, a_end] = residual_estimation(u_cheb,N,n,h,rigorous);
%     disp(['delta_0 = ',num2str(d_0)])
%     disp(['delta_tail = ',num2str(d_tail)])    
  end
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

function dy = ode_res_iter(t,y,u,gamma)

N = size(y,1);
m = (N-1)/2;

F = mapF(u, t, gamma);
buy = quadratic(u(t).',y);

k = (-m:m)';
dy = -(4*pi^2)*(k.^2).*y + 2*buy;
dy = dy*gamma + F;

end

function F = mapF(u, t, gamma)
N = (size(u,2)-1)/2;
k = (-N:N)';
du = diff(u);
u = u(t).';
u2 = quadratic(u,u);
F = du(t).' - gamma*(-(4*pi^2)*(k.^2).*u + u2);
end

% function F = mapF(u, gamma)
% N = (size(u,2)-1)/2;
% k = (-N:N);
% u2 = quad4cheb(u);
% F = diff(u) - gamma*(-(4*pi^2)*(k.^2).*u + u2);
% end

% function u2 = quad4cheb(u)
% N = (size(u,2)-1)/2;
% a = chebcoeffs(u);
% n = size(a,1);
% % get a one-sided chebyshev coefficients
% a0 = [flipud(a(2:end,:))/2;a(1,:);a(2:end,:)/2];
% a0(:,N+2:end) = fliplr(a0(:,1:N)); % symmetry-condition for Neumann boundary condition
% a = a0;  
% a2 = convtensor(a,a);
% a2 = [a2(n,:);2*a2(n+1:end,:)];
% u2 = chebfun({{[], a2}},u.domain);
% end