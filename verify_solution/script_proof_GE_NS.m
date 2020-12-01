%% preliminary
clear
% clf
% addpath('../chebfun-master/')
% addpath('../verify_defect/')
% addpath('../variational_problem/')
% addpath('../../toolbox/chebfun-master/')
addpath('../verify_defect/')
addpath('../variational_problem/')

% load ../unstable_pt.mat
% N = (size(a_unstable,1)-1)/2; % # of Fourier coefficients
N = 14; % # of Fourier coefficients
n = 13; % # of Chebyshev coefficients
% N = 15; % # of Fourier coefficients
% n = 40; % # of Chebyshev coefficients
stepsize = 0.08/2^5; % length of time step
tspan = [0,stepsize];
rigorous = 1;
% theta = pi/2;

core = 0; % core = 0 is best?

% a0 = a_unstable.';
a0 = zeros(1,2*N+1); % Initial sequence
a0(N+1)=50; a0(N)=-25; a0(N+2)=-25;
% a0(N+1)=21.7656; a0(N)=15.5256*1i; a0(N+2)=15.5256*1i; a0(N-1)=-2.0530; a0(N+3)=-2.0530; %a0(N-2)=-1i*0.2027; a0(N+4)=-1i*0.2027;
% a0(N+1)=4*21.7656; a0(N-1)=4*15.5256*1i; a0(N+3)=4*15.5256*1i; a0(N-3)=4*-2.0530; a0(N+5)=4*-2.0530; a0(N-5)=4*-1i*0.2027; a0(N+7)=4*-1i*0.2027;

if rigorous>0
  y_local = intval(zeros(1,10));
else
  y_local = zeros(1,10);
end

y = []; % Data container

%% getting approximate solution and residual bounds
% timestep = 1;
eps_all = 0;
for timestep = 1:1000
  [a, d_N, d_infty] = getting_the_solution_timestepping(N,n,tspan,a0,rigorous);% Output is one-sided Chebyshev!
  disp(['delta_N = ',num2str(d_N)])
  disp(['delta_tail = ',num2str(d_infty)])
  h = stepsize;
  
%   plot_solution(a,tspan,1), hold on, pause(0.01)
  
  % Initial error
  if rigorous>0
    eps_all = sup(eps_all + compute_eps(intval(a),intval(a0)));
    d_all = sup(d_N+d_infty);
  else
    eps_all = eps_all + compute_eps(a,a0);
    d_all = d_N+d_infty;
  end
  
  
  %% start to solve variational problem
  [C,r_minus] = solve_variational_equation(a,h,N,n,core,rigorous);
  % r_minus
  if rigorous>0
    C0 = reshape(C,n,2*core+1,2*core+1);
    C0k = reshape(sum(abs(intval(C0(1,:,:))),2),2*core+1,1);
    M_phi = max(2*(sum(abs(intval(C)))'+r_minus)-C0k);
  else
    C0 = reshape(C,n,2*core+1,2*core+1);
    C0k = reshape(sum(abs(C0(1,:,:)),2),2*core+1,1);
    M_phi = max(2*(sum(abs(C))'+r_minus)-C0k);
  end
  
  [C_backward,r_minus_backward] = solve_variational_equation_backward(a,h,N,n,core,rigorous);
  if rigorous>0
    C0_backward = reshape(C_backward,n,2*core+1,2*core+1);
    C0k_backward = reshape(sum(abs(intval(C0_backward(1,:,:))),2),2*core+1,1);
    M_psi = max(2*(sum(abs(intval(C_backward)))'+r_minus_backward)-C0k_backward);
  else
    C0_backward = reshape(C_backward,n,2*core+1,2*core+1);
    C0k_backward = reshape(sum(abs(C0_backward(1,:,:)),2),2*core+1,1);
    M_psi = max(2*(sum(abs(C_backward))'+r_minus_backward)-C0k_backward);
  end
  
  M0 = M_phi*M_psi;
  
  C0(2:end,:,:) = 2*C0(2:end,:,:);
  if rigorous>0
    phi_at_endpoint = norm(reshape(sum(intval(C0),1),2*core+1,2*core+1),1);
  else
    phi_at_endpoint = norm(reshape(sum(C0,1),2*core+1,2*core+1),1);
  end
  
  %% start to estimate the evolution operator
  if rigorous>0
    ipi = intval('pi');
%     alp_N = (core+1)^2*(2*ipi)^2*theta;
    alp_N = 0;
    a_X = a_norm(intval(a));
    a_infty = a; a_infty(:,N+1) = 0;
    a_inf_X = a_norm(intval(a_infty));
  else
%     alp_N = (core+1)^2*(2*pi)^2*theta;
    alp_N = 0;
    a_X = a_norm(a);
    a_infty = a; a_infty(:,N+1) = 0;
    a_inf_X = a_norm(a_infty);
  end
  
  if alp_N>2*a_X
    M_infinity = (1-exp(-(alp_N-2*a_X)*h))/(alp_N-2*a_X);
    M_bar_infty = (h-M_infinity)/(alp_N-2*a_X);
    kappa = 1 - 4*M0*a_inf_X^2*M_bar_infty;
  else
    M_infinity = (exp((2*a_X-alp_N)*h)-1)/(2*a_X-alp_N);
    M_bar_infty = (M_infinity-h)/(2*a_X-alp_N);
    kappa = 1 - 4*M0*a_inf_X^2*M_bar_infty;
  end
  M_infinite_sup = max(1,exp((2*a_X-alp_N)*h));
  M_infinite_at_endpoint = exp(-(alp_N-2*a_X)*h);
  
  if kappa>0
    U_matrix = [M0/kappa, 2*M0*M_infinity*a_inf_X/kappa;...
      2*M0*M_infinity*a_inf_X/kappa, M_infinite_sup+4*M0*M_infinity^2*a_inf_X^2/kappa];
    M = sup(norm(U_matrix,1));
    U_matrix_at_endpoint1 = [phi_at_endpoint*(1+2*M_psi*a_inf_X*h*U_matrix(2,1)),...
      2*phi_at_endpoint*M_psi*a_inf_X*h*U_matrix(2,2);...
      2*M_infinity*a_inf_X*U_matrix(1,1),...
      exp((2*a_X-alp_N)*h)+2*M_infinity*a_inf_X*U_matrix(1,2)];
    U_matrix_at_endpoint2 = [phi_at_endpoint*M_psi*(1+2*a_inf_X*h*U_matrix(2,1)),...
      2*phi_at_endpoint*M_psi*a_inf_X*h*U_matrix(2,2);...
      2*M_infinity*a_inf_X*U_matrix(1,1),...
      M_infinite_sup+2*M_infinity*a_inf_X*U_matrix(1,2)];
    M_at_endpoint = min(M,sup(norm(U_matrix_at_endpoint1,1)));
    Ms = min(M,sup(norm(U_matrix_at_endpoint2,1)));
  else
    error('Linearized problem is not solved (kappa<0)')
  end
  
  
  %% verify the contraction mapping: require INTLAB for using "verifynlssall"
  F = @(x) M*(eps_all+stepsize*(2*x.^2+d_all))-x;
  
  [xx,xx_cand,data] = verifynlssall(F,infsup(0,10));
  if ~any(xx>0)
    error('contraction mapping is not verified!')
  end
  while(1)
    if isempty(xx_cand)
      err = sup(xx(:,all(F(sup(xx))<0,1)));
      break
    else
      [xx,xx_cand,data] = verifynlssall(data);
    end
  end
  err_at_endpoint = intval(M_at_endpoint)*eps_all+intval(Ms)*stepsize*(2*intval(err)^2+intval(d_all));
  err_at_endpoint = min(err,sup(err_at_endpoint));
  
  %% verify global existence
  if rigorous > 0
    a0 = sum(intval(a))+sum(intval(a(2:end,:)));% initial sequence of next step
    
    
% %     phi0 = a0(N+1) + err_at_endpoint*(infsup(-1,1)+1i*infsup(-1,1)) ;
% %     tphi = a0 + err_at_endpoint*infsup(-1,1);
% %     tphi(N+1) = 0;
% %     tphi = (tphi./(phi0^2));
% % 
% %     phi2 = a0;
% %     phi2(N+1) = 0;
% %     norm_phi2 = sup(norm(phi2,1)) + err_at_endpoint;
% %     norm_tphi2 = norm_phi2/(phi0^2);
    
    phi0_center = a0(N+1);
    norm_phi0   = abs(phi0_center);
    phi0_int    = midrad( mid(phi0_center),rad(phi0_center)+err_at_endpoint);
    
    phi2        = a0;
    phi2(N+1)   = 0;
    norm_phi2   = sup(norm(phi2,1));
    
%   We have a bound 
%   | tildephi | \leq  
%         \sup_{\delta \in [0,1]}   || norm_phi2 +(1-\delta)*err_at_endpoint ||  /  ( |\phi0| - \delta*err_at_endpoint )^2
%   If the derivative with respect to delta is negative, we take delta =0,
%   and obtain  
%               | tildephi | \leq || norm_phi2 + err_at_endpoint || /|\phi0|^2
%                                                       
     if ( 2*norm_phi2 - norm_phi0 + 2*err_at_endpoint)<0
         norm_tphi2 = (norm_phi2+err_at_endpoint)/(norm_phi0)^2;
     else
         norm_tphi2 = (norm_phi2+err_at_endpoint)/(norm_phi0-err_at_endpoint)^2;
     end
    
%      Data_List(Data_List_counter,:) = [ sup(abs(phi0_center)), sup(abs(phi0_int)), sup(abs(norm_phi2/phi0_center^2)), sup(norm_tphi2)];
%      Data_List_counter = Data_List_counter+1;
    [success_GE,r] = verify_GE(phi0_int,norm_tphi2);
%     [success_GE,r] = verify_GE(phi0,tphi);
    a0 = mid(a0);
  else
    a0 = sum(a)+sum(a(2:end,:));% initial sequence of next step
    phi0 = a0(N+1);
    tphi = a0;
    tphi(N+1) = 0;
    tphi = (tphi./(phi0^2));
    [success_GE,r] = verify_GE(phi0,tphi);
  end
  
  
  %% Data
  y_local(1) = tspan(1);
  y_local(2) = tspan(2);
  y_local(3) = M0;
  y_local(4) = M_at_endpoint;
  y_local(5) = M;
  y_local(6) = a_X;
  y_local(7) = phi0_center;
  y_local(8) = eps_all;
  y_local(9) = d_all;
  y_local(10) = err;
  % y(10) = sup(min(xx));
  y = [y;y_local];
  
  if success_GE>0
    break
  end
  
  %% Update the initial error and time interval
  eps_all = err_at_endpoint;
  tspan = tspan + stepsize;% next time step
  
end

success_GE
% printresult_timestepping
% save('data_GEv1.mat','y','t')
% save data_GE_NLS.mat
