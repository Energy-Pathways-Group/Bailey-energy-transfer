% LinearPV
% function to compute Linear Ertel PV field (in 1/s^3)
% function PV = LinearPV(self,rhoBar)
%
% Linear Ertel PV  (the linear terms in Ertel PV) 1/s^3
% PV = f_0 b_z + (f0 + vorticity_z) *N^2(z)
%
% Same as above, the background PV, f0*N2, can be safely
% eliminated in constant stratification, but it remains
% important in variable stratification due to vertical
% advection.
%
% Bailey Avila        09/26/2022
% Last modified       09/26/2022

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PV = LinearPV(self,noBackground)
  % function to compute LinearPV field
 
  % get wave plus, wave minus and full field components of flow
  [u,v,w,zeta] = self.transformWaveVortexToUVWEta(self.Ap,self.Am,self.A0,self.t);
 
  dx = diff(self.x(1:2));							   % (m)
  dy = diff(self.y(1:2));							   % (m)
  dz = diff(self.z(1:2));							   % (m)
  
  % eta to rho prime
  rho_prime = zeta / (self.g/(self.rho0 * self.N2));

  % buoyancy
  b = -(self.g/self.rho0)*rho_prime;

  % full field vorticity components (note, diffZF is diffCosine, diffZG is diffSine)
  vort_z = diffX(self,v,1) - diffY(self,u,1);	% v_x - u_y	% (1/s)
     
  % buoyancy perturbation derivatives
  b_z = diffZG(self,b,1);                                 % (1/s^2)

  if exist('noBackground','var') && noBackground == 1
    PV = self.f * b_z + vort_z .* self.N2;
  else
    PV = self.f * b_z + (self.f + vort_z) .* self.N2;
  end


end
