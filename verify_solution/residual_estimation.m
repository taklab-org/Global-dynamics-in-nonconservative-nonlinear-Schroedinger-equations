function [a0, d_0, d_tail, a_end] = residual_estimation(u_cheb,N,n,h,rigorous)
% getting the residual estimation

% get a two-sided chebyshev coefficients
a = chebcoeffs(u_cheb);
% if 1
% a certain filter for app. solution
% tol = 5e-19;
% a(abs(a)<tol) = 0;
% a_end = sum(a,1);
% else
a_end = u_cheb(end);
% end

% compute a derivative of the app. sol. from the chebyshev coefficients
rescaleFactork = h/2;
du = ChebDerCoeffs(a,rigorous)/rescaleFactork;
du = [du;zeros(1,size(du,2))];
du = [du(1,:);du(2:end,:)/2];

% get a one-sided chebyshev coefficients
a0 = [flipud(a(2:end,:))/2;a(1,:);a(2:end,:)/2];
a0(:,N+2:end) = fliplr(a0(:,1:N)); % symmetry-condition for Neumann boundary condition

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
  
  ires_full(2:end,:) = 2*ires_full(2:end,:);% back to two-sided Chebyshev coefficients
  
%   d_N = mag(sum(sum(abs(ires_full(:,N+1:3*N+1)))))
%   ires_full(:,N+1:3*N+1)=0;
%   d_infty = mag(sum(sum(abs(ires_full))))
%   d_N + d_infty
  
  d_0 = sup(sum(abs(ires_full(:,2*N+1))));
  ires_full(:,2*N+1)=0;
  d_tail = sup(sum(sum(abs(ires_full))));
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
  res_full = du_full + 1i*((4*pi^2)*(k_full.^2).*a_full-a2full);
  
%   res(2:end,:) = 2*res(2:end,:);% back to chebyshev coefficients
  res_full(2:end,:) = 2*res_full(2:end,:);% back to chebyshev coefficients
%   disp(sum(sum(abs(res),2)))
%   disp(sum(sum(abs(res_full),2)))
%   disp(sum(sum(abs(res_full),2))-sum(sum(abs(res),2)))
%   d_N = sum(sum(abs(res_full(:,N+1:3*N+1))));
%   res_full(:,N+1:3*N+1)=0;
%   d_infty = sum(sum(abs(res_full)));
  d_0 = sup(sum(abs(res_full(:,2*N+1))));
  res_full(:,2*N+1)=0;
  d_tail = sup(sum(sum(abs(res_full))));
end
