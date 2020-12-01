% Exectime almost 921.033635 sec. by destop PC (Core i9 10900K, Nov.20, 2020)
%% preliminary
clear
clc
addpath('../verify_defect/')
addpath('../variational_problem/')

load end_points_proof_NLS_pt2_conj
% load end_points_proof_NLS_pt2_v5
N = (size(P_at_1,1)-1); % # of Fourier projection
n = 7; % # of Chebyshev coefficients % 2^4:9, 2^5:8, 2^6:7
stepsize = 0.08/2^7; % initial length of time step
tspan = [0,stepsize];
core = 0; % core = 0 is best

a0 = mid([flipud(P_at_1(2:end));P_at_1].');
err_at_endpoint_old = sup([error_at_1; error_at_1]);

max_iteration = 1e4;

y_local = intval(zeros(1,10));
y = []; % Data container

initial_step = true;

%% getting approximate solution and residual bounds
for timestep = 1:max_iteration
  [a, d_0, d_tail] = getting_the_solution_timestepping_alt(N,n,tspan,a0);% Output is one-sided Chebyshev!
  disp(['delta_0 = ',num2str(d_0)])
  disp(['delta_tail = ',num2str(d_tail)])
  h = stepsize;
  
  %     plot_solution(a,tspan,1), hold on, pause(0.01)
  
  % Initial error
  eps_all = sup(err_at_endpoint_old + compute_eps_alt(intval(a),intval(a0)));
  d_all = sup([d_0;d_tail]);
  
  %% verify local existence
  [err,M_at_endpoint,Ms,M0,M,a_X,~] = verify_local_existence(eps_all,d_all,a,h,N,n,core);
  if any(isnan(err))
    break
  end
  
  %% estimate err_at_endpoint and adjust time step
  err_at_endpoint = intval(M_at_endpoint)*eps_all+intval(Ms)*stepsize*((2*sum(intval(err)))^2*ones(2,1)+intval(d_all));
  err_at_endpoint = min(err,sup(err_at_endpoint));
  
  % adjust stepsize corresponding to increase ratio of err_at_endpoint
  if any(err_at_endpoint./err_at_endpoint_old>1.05) %&& ~initial_step
    stepsize = stepsize/(max(err_at_endpoint./err_at_endpoint_old))^2;
    tspan(2) = tspan(1) + stepsize;% next time step
    disp('adjust timestep (smaller)')
    continue
  elseif all(err_at_endpoint./err_at_endpoint_old <= 1.001)
    stepsize = stepsize*1.1;
    tspan(2) = tspan(1) + stepsize;% next time step
    disp('adjust timestep (larger)')
    continue
  end
  
  
  %% verify global existence
  a0 = sum(intval(a))+sum(intval(a(2:end,:)));% initial sequence of next step
  
  
  phi0_center = a0(N+1);
  phi0_int    = midrad( mid(phi0_center),0);
  norm_phi0   = abs(phi0_int);
  
  
  phi2        = a0;
  phi2(N+1)   = 0;
  norm_phi2   = sup( err_at_endpoint(1) + rad(phi0_center) ...
    +  err_at_endpoint(2)+norm(phi2,1) );
  norm_tphi2   = norm_phi2 / norm_phi0^2;
  
  disp([ norm_phi0 abs(norm_tphi2)])
  phi0_int
  tspan
  
  [success_GE,r] = verify_GE(phi0_int,norm_tphi2);
  %     [success_GE,r] = verify_GE(phi0,tphi);
  a0 = mid(a0);
  initial_step = false;
  
  %% save data
  y_local(1) = tspan(1);
  y_local(2) = tspan(2);
  y_local(3) = M0;
  y_local(4) = norm(M_at_endpoint,1);
  y_local(5) = norm(M,1);
  y_local(6) = a_X;
  y_local(7) = phi0_center;
  y_local(8) = norm(err_at_endpoint,1);
  y_local(9) = norm(d_all,1);
  y_local(10) = norm(err,1);
  y = [y;y_local];
  
  if success_GE>0
    return
  end
  
  %% Update the initial error and time interval: old means the previous step
  err_at_endpoint_old = err_at_endpoint;
  tspan(1) = tspan(2);
  tspan(2) = tspan(2) + stepsize;% next time step
end

success_GE
% printresult_timestepping
% save('data_GEv1.mat','y','t')
% save data_NLS_from_P_at_minus_1.mat
