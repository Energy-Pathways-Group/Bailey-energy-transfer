% APV
% function to compute APV field (APV in 1/s) per Early et al (in process)
% function PV = APV(self)
%
% APV = (\nabla \times u) \cdot \hat{z} + f_0 \frac{\partial \eta}{\partial z} - (\nabla \times u) \cdot \nabla \eta
%
% 
% Bailey Avila        12/05/2022
% Last modified       12/05/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PV = AvailPV_early(self)
  % function to compute Available PV field per Early et al paper 
 
  % get wave plus, wave minus and full field components of flow
  [u,v,w,eta] = self.transformWaveVortexToUVWEta(self.Ap,self.Am,self.A0,self.t);

  % full field vorticity components (note, diffZF is diffCosine, diffZG is diffSine)
  vort_x = diffY(self,w,1) - diffZF(self,v,1);	% w_y - v_z	% (1/s)
  vort_y = diffZF(self,u,1) - diffX(self,w,1);	% u_z - w_x	% (1/s)
  vort_z = diffX(self,v,1) - diffY(self,u,1);	% v_x - u_y	% (1/s)
     
  % buoyancy perturbation derivatives
  eta_x = diffX(self,eta,1);                                  % (1/s^2)
  eta_y = diffY(self,eta,1);                                  % (1/s^2)
  eta_z = diffZG(self,eta,1);                                 % (1/s^2)

  % now sum components to get PV
  PV = (vort_z) - self.f*eta_z - (vort_x.*eta_x) - (vort_y.*eta_y) - (vort_z.*eta_z);

end
