function [success,r] = verify_GE(phi0,tphi)
success = 0;
if real(1i*phi0) < 0
  ipi = sup(intval('pi'));
  norm_phi0 = sup(abs(phi0));
  norm_tphi = sup(norm(tphi,1));
  
  alpha = ipi*norm_phi0/2;
  beta  = (1/norm_tphi - alpha);

  discriminant =  beta^2 - 2 * alpha^2;
  if beta < 0 || discriminant < 0
      r = NaN;
      disp('Endpoint is too large to verify global existence.')
      return
  else 
      delta = sqrt(discriminant);
  end
  
  r_min = (beta - delta)/(alpha^2);
  r_max = (beta + delta)/(alpha^2);
  r_mid = r_min*1.05;
  
  if norm_tphi * exp( alpha * r_mid) < r_mid
      r = r_mid;
      success = 1; 
      return
  end
      
  
% %   Radii_search_space = infsup(0,10)*ones(2,1);

% %   r1_min = (pi/2)*norm_tphi^2*norm_phi0;
% %   r2_min = norm_tphi;
% %   
% %   r1_max = (norm_phi0*norm_tphi^2)^-1 - pi^2/8;
% %   r2_max = (norm_tphi^-1 - norm_phi0*pi/2)/( (norm_phi0*pi/2)^2*3/2 );
% %   
% %   delta  = (log(r2_max/norm_tphi)-r2_min*norm_phi0*pi/2)/(2*r2_min^2*norm_phi0^2);
% %   
% %   r1_max2 = (delta - pi^2/8)/(pi^3/48);
  
  r_max = min([r_max,10]); 
  
  disp(['Search space:  r1 = [', num2str(r_min) ' , ' num2str(r_max) ' ]'])
  
%   if ~( (r1_min < r1_max) && (r2_min < r2_max) )
%       r = NaN;
%       disp('Endpoint is too large to verify global existence.')
%       return
%   end
  
  r_search_space = infsup(r_min,r_max);
  Radii_search_space = [ r_search_space ];
  
% % % % %   if norm_phi0 > .875
% % % % %       r = -intval(ones(2,1));
% % % % %       return
% % % % %   end
  G = @(r) norm_tphi * exp( alpha * r(1,:)) - r(1,:);


%   G = @(r) [ ( r(2,:).^2 *  norm_phi0 ./ r(1,:) ) .* (exp(r(1,:)*ipi/2 ) -1) - r(1,:);...
%      norm_tphi * exp(...
%      ipi/2 .* r(2,:).*norm_phi0...
%      + (2*norm_phi0.^2 .* r(2,:).^2 ./ r(1,:))...
%      .*( (exp(r(1,:).*ipi/2) -1)./r(1,:) - ipi/2  ) ...
%      ) - r(2,:)];
  [r,r_cand,data] = verifynlssall(G,Radii_search_space);
  if ~any(r>0)
    disp('global existence is not verified.')
  end
  while(1)
    if isempty(r_cand)
      if isempty(r)
        break
      end
      r = sup(r(:,all(G(sup(r))<0,1)));
      success = 1;
      break
    else
      [r,r_cand,data] = verifynlssall(data);
    end
  end
else
  r = nan;
  disp('Re(1i*phi_0) is non-negative.')
end
