function M = getting_M_t1s(C0,C0_backward,r_minus,r_minus_backward,h,n,core)
% Back to two-sided chebyshev
C0(2:end,:,:) = 2*C0(2:end,:,:);
C0_backward(2:end,:,:) = 2*C0_backward(2:end,:,:);

size_of_matrix = size(C0,3);
% n = 13; % number of chebyshev polynomial

%%%%%%%%%%% Here should be computed one time %%%%%%%%%%%%%%%
% Creat a mesh based on Chebyshev points of higher order
m = 2^5; % divide number
theta = linspace(intval('pi'),0,m*n+1); % angles of chebyshev points divided by m
t = 0.5*h*(1+cos(theta)); % values of mesh points
mesh = hull(t(1:end-1),t(2:end)); % construct the mesh
% mesh_size = length(mesh);
% y = [];
% for i = 1:n % chebyshev nodes
% cheb_val = cos(i*acos(-(h-2*mesh)/h)); % values of chebyshev polynomials at mesh domain
% y = [y;[abs(C0(i,1,1)), max(abs(C0(i,1,1)*cheb_val))]];
% end
% [max(abs(C0(1:n,1,1).*cheb_val),[],2),y]

% Chebyshev values on the mesh
i = (1:n)';
s_cheb_val = (cos(i*acos(-(h-2*mesh)/h)));

% % Create indices for mesh
% aa = 1:mesh_size; e = ones(1,mesh_size);A = [kron(aa,e);kron(e,aa)];
% dom_ind = A(:,A(1,:)>=A(2,:)); % dom_ind(1,:) -> t_index, dom_ind(2,:) -> s_index
% 
% t_cheb_val = cheb_val(:,dom_ind(1,:));
% s_cheb_val = cheb_val(:,dom_ind(2,:));

%%%%%%%%%%%%%%%%%%%%%%%%%%

num_nodes = length(s_cheb_val);

% Phi values on each mesh
% CC0 = (reshape(C0,n,size_of_matrix^2).')*t_cheb_val + r_minus*infsup(-1,1);
% Phi_value_on_mesh = reshape(CC0.',num_nodes,size_of_matrix,size_of_matrix);
Phi_value_on_mesh = repmat(sum(C0,1) + r_minus*infsup(-1,1),num_nodes,1,1);

% Psi values on each mesh
CC0_bakckward = (reshape(C0_backward,n,size_of_matrix^2).')*s_cheb_val + r_minus_backward*infsup(-1,1);
Psi_value_on_mesh = reshape(CC0_bakckward.',num_nodes,size_of_matrix,size_of_matrix);

Phi_Psi_value_on_mesh = intval(zeros(num_nodes,size_of_matrix,size_of_matrix));
% tic
for i=1:2*core+1
  for j=1:2*core+1
     Phi_Psi_value_on_mesh(:,i,j) = sum(reshape(Phi_value_on_mesh(:,i,:),num_nodes,size_of_matrix).*Psi_value_on_mesh(:,:,j),2);
  end
end
% toc
% M0 = M_phi*M_psi

M = max(max(sum(abs(Phi_Psi_value_on_mesh),2),[],3));
