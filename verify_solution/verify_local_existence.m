function [err,W_at_endpoint,W_J,Wm,W_h,ba_X,kappa] = verify_local_existence(eps_all,d_all,a,h,N,n,core,rigorous)
arguments
  eps_all; d_all; a; h; N; n; core;
  rigorous = true;
end
%% start solving variational problem
  [C,r_minus] = solve_variational_equation(a,h,N,n,core,rigorous);
  % r_minus
  if rigorous
    C0 = reshape(C,n,2*core+1,2*core+1);
    C0k = reshape(sum(abs(intval(C0(1,:,:))),2),2*core+1,1);
    M_phi = max(2*(sum(abs(intval(C)))'+r_minus)-C0k);
  else
    C0 = reshape(C,n,2*core+1,2*core+1);
    C0k = reshape(sum(abs(C0(1,:,:)),2),2*core+1,1);
    M_phi = max(2*(sum(abs(C))'+r_minus)-C0k);
  end
  
  [C_backward,r_minus_backward] = solve_variational_equation_backward(a,h,N,n,core,rigorous);
  if rigorous
    C0_backward = reshape(C_backward,n,2*core+1,2*core+1);
    C0k_backward = reshape(sum(abs(intval(C0_backward(1,:,:))),2),2*core+1,1);
    M_psi = max(2*(sum(abs(intval(C_backward)))'+r_minus_backward)-C0k_backward);
  else
    C0_backward = reshape(C_backward,n,2*core+1,2*core+1);
    C0k_backward = reshape(sum(abs(C0_backward(1,:,:)),2),2*core+1,1);
    M_psi = max(2*(sum(abs(C_backward))'+r_minus_backward)-C0k_backward);
  end
% % A technique to estimate a uniform bound using interval arithmetic
  if rigorous
    Wm = getting_M0(C0,C0_backward,r_minus,r_minus_backward,h,n,core);%sup(Wm)
    Wmt = getting_M_t1s(C0,C0_backward,r_minus,r_minus_backward,h,n,core);
    Wm = min(M_phi*M_psi,Wm);
    Wmt = min(M_phi*M_psi,Wmt);
  else
    Wm = sup(getting_M0(C0,C0_backward,r_minus,r_minus_backward,h,n,core));
    Wm = min(M_phi*M_psi,Wm);
    Wmt = sup(getting_M_t1s(C0,C0_backward,r_minus,r_minus_backward,h,n,core));
    Wmt = min(M_phi*M_psi,Wmt);
  end
  
  C0(2:end,:,:) = 2*C0(2:end,:,:);
  if rigorous
    phi_at_endpoint = norm(reshape(sum(intval(C0),1),2*core+1,2*core+1),1);
  else
    phi_at_endpoint = norm(reshape(sum(C0,1),2*core+1,2*core+1),1);
  end
  
  %% start estimating the evolution operator
  if rigorous
    mu_m = 0;
    ba_X = a_norm(intval(a));
    a_infty = a; a_infty(:,N+1) = 0;
    ba_inf_X = a_norm(intval(a_infty));
    ba_inf_X_dual = max(sum(abs(intval(a_infty)),1));
%     ba_inf_X_dual = a_norm(intval(a_infty)); % previous estimate
  else
    mu_m = 0;
    ba_X = a_norm(a);
    a_infty = a; a_infty(:,N+1) = 0;
    ba_inf_X = a_norm(a_infty);
    ba_inf_X_dual = max(sum(abs(a_infty),1));
  end
  
  if mu_m>2*ba_X
    W_infinity = (1-exp(-(mu_m-2*ba_X)*h))/(mu_m-2*ba_X);
    W_bar_infty = (h-W_infinity)/(mu_m-2*ba_X);
    kappa = 1 - 4*Wm*ba_inf_X_dual*ba_inf_X*W_bar_infty;
  else
    W_infinity = (exp((2*ba_X-mu_m)*h)-1)/(2*ba_X-mu_m);
    W_bar_infty = (W_infinity-h)/(2*ba_X-mu_m);
    kappa = 1 - 4*Wm*ba_inf_X_dual*ba_inf_X*W_bar_infty;
  end
  W_infinite_sup = max(1,exp((2*ba_X-mu_m)*h));
%   W_infinite_at_endpoint = exp(-(alp_N-2*a_X)*h);
  
  if kappa>0
    U_matrix = [Wm, 2*Wm*W_infinity*ba_inf_X_dual;...
      2*Wm*W_infinity*ba_inf_X, W_infinite_sup]/kappa;
    U_matrix_at_endpoint1 = [phi_at_endpoint+2*Wmt*ba_inf_X_dual*h*U_matrix(2,1),...
      2*Wmt*ba_inf_X_dual*h*U_matrix(2,2);...
      2*W_infinity*ba_inf_X*U_matrix(1,1),...
      exp((2*ba_X-mu_m)*h)+2*W_infinity*ba_inf_X*U_matrix(1,2)];
    U_matrix_J = [Wmt*(1+2*ba_inf_X_dual*h*U_matrix(2,1)),...
      2*Wmt*ba_inf_X_dual*h*U_matrix(2,2);...
      2*W_infinity*ba_inf_X*U_matrix(1,1),...
      W_infinite_sup+2*W_infinity*ba_inf_X*U_matrix(1,2)];
  else
    error('Linearized problem is not solved (kappa<0)')
  end
  
  
  %% verify the contraction mapping: require INTLAB for using "verifynlssall"
  if length(eps_all)>1
    U_matrix_sup = sup(U_matrix);
    F = @(x) U_matrix_sup*[eps_all(1)+h*(2*(x(1,:).^2+x(2,:).^2)+d_all(1));eps_all(2)+h*(2*(x(1,:)+x(2,:)).^2+d_all(2))]-[x(1,:);x(2,:)];
    finite_infinite = 2;
    W_at_endpoint = U_matrix_at_endpoint1;
    W_J = U_matrix_J;
    W_h = U_matrix;
  else
    W_h = sup(norm(U_matrix,1));
    F = @(x) W_h*(eps_all+h*(2*x.^2+d_all))-x;
    finite_infinite = 1;
    W_at_endpoint = min(W_h,sup(norm(U_matrix_at_endpoint1,1)));
    W_J = min(W_h,sup(norm(U_matrix_J,1)));
  end

  
  [xx,xx_cand,data] = verifynlssall(F,infsup(0,1e2)*ones(finite_infinite,1));
  if ~any(xx>0)
    disp('contraction mapping is not verified!')
    err = NaN;
    return
  end
  while(1)
    if isempty(xx_cand)
      err = sup(xx(:,all(F(sup(xx))<0,1)));
      break
    else
      [xx,xx_cand,data] = verifynlssall(data);
    end
  end