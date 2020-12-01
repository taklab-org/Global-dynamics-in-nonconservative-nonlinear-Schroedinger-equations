clear
close all
load ../verify_solution/data_NLS_from_P_at_1.mat
% load ../verify_solution/data_NLS_from_P_at_minus_1.mat % gives err_bound_5_2 & stepsize_5.2
% load ../verify_solution/data_NLS_pt2.mat % gives err_bound_5_3 & stepsize_5_3
% load ../verify_solution/data_NLS_pt1_conj_from_P_at_1.mat
% load ../verify_solution/data_NLS_pt1_conj_from_P_at_minus_1.mat
% load ../verify_solution/data_NLS_pt2_conj.mat


LW = 'linewidth'; lw = 1.6;
FS = 'FontSize';% fs = 20;

plot(y(:,2).mid- y(:,1).mid, LW, lw)
set(gca,'FontSize',15)
xlabel('Number of time steps',FS,20)
ylabel('$$h_i$$','interpreter', 'latex',FS,25)
set(gcf, 'position', [380 320 650 350]);
title('Adaptive step size',FS,20)

figure 
semilogy(mid(y(:,2)),mid(y(:,8)),LW,lw), hold on
% semilogy(mid(y(:,2)),mid(y(:,10)),LW,lw), hold on
xlabel('$$t$$','interpreter', 'latex',FS,25)
ylabel('$$\varrho_0+\varrho_\infty$$','interpreter', 'latex',FS,25)
set(gcf, 'position', [380 320 650 350]);
title('Rigorous error bounds',FS,20)