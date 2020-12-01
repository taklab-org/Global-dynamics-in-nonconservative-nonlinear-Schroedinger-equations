%%
clear
close all
% addpath('../toolbox/chebfun-master/')
% addpath('../verify_solution/')
% addpath('../verify_defect/')
% pic = 'P_at_1';
pic = 'P_at_minus_1';
% pic = 'pt2';

switch pic
  case 'P_at_1' % pt1_P_at_1
    %% CO_pt1_from_P_at_1
    addpath('../Manifolds/')
    load unstable_manifold1
    f1 = figure;
    plot_manifold_v2(p,M,N,lambda,1);
    rmpath('../Manifolds/')
    load end_points_proof2_NLS.mat
    N = (size(P_at_1,1)-1); % # of Fourier projection
    a0 = ([flipud(P_at_1(2:N+1));P_at_1(1:N+1)].');
    prooftime = 0.05;
    endtime = 0.2;

    plot_integration(f1,a0,N,endtime)
    figure(f1)
    subplot(2,1,1)
    view([90, 90])
    colorbar
    title('$$\mathrm{Re}(u_a)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[prooftime,prooftime],[50,50],'w-','linewidth',2)
    subplot(2,1,2)
    view([90, 90])
    colorbar
    title('$$\mathrm{Im}(u_a)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[prooftime,prooftime],[50,50],'k-','linewidth',2)
    f1.Position = [380 320 1097 368*2];
    colormap jet

    %% CO_pt1_conj_from_P_at_1
    addpath('../Manifolds/')
    load unstable_manifold1_conj
    f2 = figure;
    plot_manifold_v2(p,M,N,lambda,1,true);
    rmpath('../Manifolds/')
    load end_points_proof_NLS_pt1_conj
    N = (size(P_at_1,1)-1); % # of Fourier projection
    a0 = [flipud(P_at_1(2:N+1));P_at_1(1:N+1)].';
    prooftime = 0.05;
    endtime = 0.2;

    plot_integration(f2,a0,N,endtime,true)
    figure(f2)
    subplot(2,1,1)
    view([90, 90])
    colorbar
    title('$$\mathrm{Re}(u_b)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[-prooftime,-prooftime],[50,50],'w-','linewidth',2)
    subplot(2,1,2)
    view([90, 90])
    colorbar
    title('$$\mathrm{Im}(u_b)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[-prooftime,-prooftime],[50,50],'k-','linewidth',2)
    f2.Position = [380 320 1097 368*2];
    colormap jet

  case 'P_at_minus_1' % pt1_P_at_minus_1
    %% CO_pt1_from_P_at_minus_1
    addpath('../Manifolds/')
    load unstable_manifold1
    f1 = figure;
    plot_manifold_v2(p,M,N,lambda,-1);
%     [f1,f2] = plot_manifold_v2(p,M,N,lambda,-1);
    rmpath('../Manifolds/')
    load end_points_proof2_NLS.mat
    N = (size(P_at_minus_1,1)-1); % # of Fourier projection
    a0 = ([flipud(P_at_minus_1(2:N+1));P_at_minus_1(1:N+1)].');
    prooftime = 0.908344227036939;
    endtime = 1;

    
    plot_integration(f1,a0,N,endtime)
    figure(f1)
    subplot(2,1,1)
    view([90, 90])
    colorbar
    title('$$\mathrm{Re}(u_a)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[prooftime,prooftime],[50,50],'w-','linewidth',2)
    subplot(2,1,2)
    view([90, 90])
    colorbar
    title('$$\mathrm{Im}(u_a)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[prooftime,prooftime],[50,50],'k-','linewidth',2)
    f1.Position = [380 320 1097 368*2];
    colormap jet

    
    %% CO_pt1_conj_from_P_at_minus
    addpath('../Manifolds/')
    load unstable_manifold1_conj
    f2 = figure;
    plot_manifold_v2(p,M,N,lambda,-1,true);
    rmpath('../Manifolds/')
    load end_points_proof_NLS_pt1_conj
    N = (size(P_at_minus_1,1)-1); % # of Fourier projection
    a0 = ([flipud(P_at_minus_1(2:N+1));P_at_minus_1(1:N+1)].');
    prooftime = 0.908437416438651;
    endtime = 1;

    plot_integration(f2,a0,N,endtime,true)
    figure(f2)
    subplot(2,1,1)
    view([90, 90])
    colorbar
    title('$$\mathrm{Re}(u_b)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[-prooftime,-prooftime],[50,50],'w-','linewidth',2)
    subplot(2,1,2)
    view([90, 90])
    colorbar
    title('$$\mathrm{Im}(u_b)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[-prooftime,-prooftime],[50,50],'k-','linewidth',2)
    f2.Position = [380 320 1097 368*2];
    colormap jet

  case 'pt2' % pt2_P_at_1
    %% CO_pt2
    load end_points_proof_NLS_pt2
    addpath('../Manifolds/')
    load unstable_manifold2_v2
    f1 = figure;
    plot_manifold_v2(p,M,N,lambda,1);
    N = (size(P_at_1,1)-1); % # of Fourier projection
    a0 = [flipud(P_at_1(2:N+1));P_at_1(1:N+1)].';
    rmpath('../Manifolds/')
    prooftime = 0.042004641313433;
    endtime = 0.1;
    
    plot_integration(f1,a0,N,endtime)
    figure(f1)
    subplot(2,1,1)
    view([90, 90])
    colorbar
    title('$$\mathrm{Re}(u_a)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[prooftime,prooftime],[50,50],'w-','linewidth',2)
    subplot(2,1,2)
    view([90, 90])
    colorbar
    title('$$\mathrm{Im}(u_a)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[prooftime,prooftime],[50,50],'k-','linewidth',2)
    f1.Position = [380 320 1097 368*2];
    colormap jet
%     f1.Position = [380 320 1097 368];
%     colormap jet
    
    
    %% CO_pt2_conj
    load end_points_proof_NLS_pt2_conj
    addpath('../Manifolds/')
    load unstable_manifold2_conj
    f2 = figure;
    plot_manifold_v2(p,M,N,lambda,1,true);
    N = (size(P_at_1,1)-1); % # of Fourier projection
    a0 = [flipud(P_at_1(2:N+1));P_at_1(1:N+1)].';
    % [f1,f2] = plot_manifold_v2(p,M,N,lambda,-1);
    % N = (size(P_at_minus_1,1)-1); % # of Fourier projection
    % a0 = [flipud(P_at_minus_1(2:end));P_at_minus_1].';
    rmpath('../Manifolds/')
    prooftime = 0.042002666962546;
    endtime = 0.1;
        
    plot_integration(f2,a0,N,endtime,true)
    figure(f2)
    subplot(2,1,1)
    view([90, 90])
    colorbar
    title('$$\mathrm{Re}(u_b)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[-prooftime,-prooftime],[50,50],'w-','linewidth',2)
    subplot(2,1,2)
    view([90, 90])
    colorbar
    title('$$\mathrm{Im}(u_b)$$', 'Interpreter', 'latex', 'FontSize', 30)
    plot3([-2,2],[-prooftime,-prooftime],[50,50],'k-','linewidth',2)
    f2.Position = [380 320 1097 368*2];
    colormap jet
    
end


%%

function plot_integration(f,a0,N,endtime,conj)
arguments
  f;a0;N;endtime;conj=false;
end

n = 20; % # of Chebyshev coefficients
stepsize = 0.08/2^3; % length of time step
tspan = [0,stepsize];
gamma = 1i;

phi0 = a0(N+1);

y=[];
t=[];


%% getting approximate solution
for timestep = 1:1e4
chebfunpref.setDefaults('factory');
chebfunpref.setDefaults('fixedLength',n);
  %
  opts = odeset('abstol',1e-18,'reltol',2.22045e-14);
  u_cheb = chebfun.ode45(@(t,y) ode_func(y,gamma),tspan,a0,opts);
  ta = chebcoeffs(u_cheb);
  ta = [ta(1,:);0.5*ta(2:end,:)];
  
  plot_glue_manifold(ta,tspan,f,conj,200)%, pause(0.01)
  
  a0 = u_cheb(end);% initial sequence of next step
  tspan = tspan + stepsize;% next time step
  if tspan(2)>endtime
    break
  end
chebfunpref.setDefaults('factory');
end
end

function dy = ode_func(y,gamma)

N = size(y,1);
fy = quadratic(y,y);
m = (N-1)/2;
k = (-m:m)';
dy = -(4*pi^2)*(k.^2).*y+fy;
dy = dy*gamma;

end