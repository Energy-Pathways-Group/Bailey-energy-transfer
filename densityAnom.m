% Single density anomaly simulation
% Using Jeffrey's new waveVortexModel (specifically WVTransform), initialize wtih 
% one density anomaly.
% 
%
% Bailey Avila        01/11/2023
% Last modified       01/11/2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Set up directories and flags
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% setup!
clear all                                          % clear everything in workspace

addpath(genpath('~/matlab'));
addpath(genpath('/Users/baileyjavila/Documents/SMAST/IWVM/GLOceanKit-integrator-tracer-floats-refactor'));
addpath(genpath('/Users/baileyjavila/Documents/SMAST/IWVM/GLNumericalModelingKit-master'));

%basedir = '/home/bavila/IWVM/';
%basedir = '/usr3/projects/IWVM/';

%runroot = 'GM_500_01r07IWr02IW';
%runroot = 'GM_500_01r07IWr01_BPIW';

runDIR = ['/Users/baileyjavila/Library/CloudStorage/OneDrive-UniversityofMassachusettsDartmouth/Coastal PO/'];
%plotDIR = [basedir,'model_raw/',runroot,'/50_thresh/'];
%plotDIR = [basedir,'model_processed/',runroot,'/'];         % assign directory to save plots

% make model_processed directory if it does not already exist
%if ~exist(plotDIR)~=0
%  mkdir(plotDIR);
%end

movieflag = 0;
printflag = 0;
energyExchange = 1; % toggle if we allow energy to go between IWs and VM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 256;

Lx = 42959;
Ly = Lx/8;
Lz = 500;

Nx = N;
Ny = N/8;
Nz = (N/2)+1; % 2^n + 1 grid points, to match the Winters model, but 2^n ok too.

f=8.0707E-5;

latitude = asind(f/2/7.2921E-5);
N0 = 5.2e-3; % Choose your stratification 7.6001e-04
g = 9.81; % (m/s^2)
rho_0 = 1025; % (kg/m^3)

% create x,y,z arrays for later use (in case not running Early wavemodel
x = [0:Lx/Nx:Lx-Lx/Nx]';
y = [0:Ly/Ny:Ly-Ly/Ny]';
z = [0:Lz/(Nz-1):Lz]' - 500;
[Y,X,Z] = meshgrid(y,x,z);

% compute inertial period
inert_per = 2*pi/f;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0,latitude=latitude);

noBackground = 1;
outputVar = WVVariableAnnotation('LinearPV',{'x','y','z'},'1/s^3', 'linear potential vorticity');
f1 = @(wvt) LinearPV(wvt,noBackground);
wvt.addOperation(WVOperation('LinearPV',outputVar,f1));

outputVar2 = WVVariableAnnotation('AvailPV_early',{'x','y','z'},'1/s', 'available potential vorticity early');
f2 = @(wvt) AvailPV_early(wvt);
wvt.addOperation(WVOperation('AvailPV_early',outputVar2,f2));

% initialize arrays (zeros to have appropriate size)
[u, v, s1] = meshgrid(zeros(size(y)), zeros(size(x)), zeros(size(z)));

% anomaly specs in x and z
anom_x = 5*(Lx/Nx);
anom_z = 5*(Lz/(Nz-1));

% Calculate density from stratification, or else just specify rhobar directly.
z = linspace(-Lz,0,Nz);
rhobar = -rho_0/g * cumtrapz(-z, N0^2*ones(size(z)));	% N0^2 = -g/rho_0*d/dz(rhobar), given N0^2, solve for rhobar
rhobar = rhobar - rhobar(end);					% shift rhobar to make positive for all z, don't change d/dz(rhobar)
% this makes rhobar in this workspace agree w/ wavemodel rhobar

% for arbitrary stratification, have to specify a function handle.
rhobar = rhobar + rho_0;					% make rhobar include rho_0

s1 = 15*(rhobar(3)-rhobar(2)).*exp( -(-Z - (Lz/2)).^2 / (2*anom_z^2) ) ...
    .*  exp(-((X-(Lx/2)).^2 + (Y-(Ly/2)).^2) / (2*anom_x^2) );   % set single density anomaly

% go from s1 to eta
eta = s1 * g/(rho_0 * N0.^2);

% initialize model
wvt.initWithUVEta(u,v,eta);

%keyboard


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nonlinear operator
%nlOp = NonlinearBoussinesqWithReducedInteractionMasks(wvt);

% do/do not allow for energy exchange between internal waves and vortical model
if energyExchange==0
  nlOp = NonlinearBoussinesqWithReducedInteractionMasks(wvt);
  nlOp.freezeEnergyOfConstituents(WVFlowConstituents.geostrophic); % does not allow energy exchange 
else
  %nlOp = NonlinearBoussinesqWithReducedInteractionMasks(wvt);
  nlOp = BoussinesqConstantN(wvt,shouldUseSpectralVanishingViscosity=1,shouldAntialias=1,uv_damp=wvt.uMax,w_damp=wvt.wMax);
end

% initialize the integrator with the model
model = WVModel(wvt,nonlinearFlux=nlOp); % nonlinear run
%model = WVModel(wvt); % linear run

% set initial positions for a bunch of floats
nTrajectories = 100;%101; % number of floats
xFloat = Lx/2*ones(1,nTrajectories);
yFloat = Ly/2*ones(1,nTrajectories);
zFloat = linspace(-Lz,0,nTrajectories);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Np = sqrt(nTrajectories);            % N should be an even number
nFloats = nTrajectories;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make narrow 'streak' of drifters at middle of x-domain, keeping careful track of indices of drifter locations
% Note: Drifters must be on x,y grid for rho interpolation to work correctly (??, but perhaps not? - see below)

xind = Nx/2 + [-Np/2+1 : Np/2];                    % put drifters along center of domain @ even gridpoint separations
xf = x(xind);

% distribute N+1 drifters evenly along y-domain
yf = linspace(0, Ly, Np+1);                       % add extra row, but then ...
yf = yf(1:Np);                                       % drop last one as duplicate on periodic domain
% adjust y-position to be on nearest gridpoint
for nn=1:length(yf)
  [junk,yind(nn)] = min(abs(y-yf(nn)));
  yf(nn) = y(yind(nn));
end

% center of dye streak at middle of z-domain - note, Jeffrey assumes z=0 at surface, so zf should be negative here
zf = -Lz/2;

% create full array of positions
[xFloat,yFloat,zFloat] = ndgrid(xf,yf,zf);

xFloat = xFloat(:);
yFloat = yFloat(:);
zFloat = zFloat(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now place drifters on an isopycnal.
% last argument can be 'int' or 'both' to use internal modes only, or internal + external
% Note: z_float should be inherently negative here, to be consistent with Jeffrey's convention
% Good for GLOceanKit_012018; should agree w/ internal modes only conditional when initializing u,v,w velocity fields
disp(' Placing particles on isopycnals using internal modes only ...')
%[z_isopycnal, rho_isopycnal] = wavemodel.PlaceParticlesOnIsopycnal(xFloat,yFloat,zFloat, ...
%  'interpolationMethod',interpolationMethod, 'tolerance',1e-8, 'maxIterations',200, ...
%  'useModes','internalOnly', 'shouldShowDiagnostics',true);
model.setFloatPositions(xFloat,yFloat,zFloat,'AvailPV_early','LinearPV');

% make x_float, y_float, z_float now contain float positions on isopycnals.
%zFloat = z_isopycnal;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At this point have checked everything regarding floats, so do the final step of converting from Jeffrey's z-coordinate
% of -Lz:0 to Winters' coordinate of 0:Lz.
%keyboard
% Shift vertical position to flow_solve coordinates (zero at bottom)
zFloat = zFloat + Lz;

%model.setFloatPositions(xFloat,yFloat,zFloat,'rho_total');

%modelDT = period/100;
% Set up the integrator
%nT = model.setupIntegrator(timeStepConstraint="min", outputInterval=period/15,finalTime=3*period,cfl=0.25);
%nT = model.setupIntegrator(outputInterval=period/10,finalTime=3*period,deltaT=modelDT);
%nT = model.setupIntegrator(timeStepConstraint="oscillatory", outputInterval=period/10,finalTime=3*period);
%nT = model.setupIntegrator(timeStepConstraint="min", outputInterval=period/15,finalTime=3*period);

%nT = model.setupIntegrator(timeStepConstraint="oscillatory", outputInterval=period/10,finalTime=3*period);
%nT = model.setupIntegrator(timeStepConstraint="oscillatory", outputInterval=150,deltaT=12,finalTime=864000);
nT = model.setupIntegrator(timeStepConstraint="oscillatory", outputInterval=180,deltaT=12,finalTime=864000/10);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Track floats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% track floats
model.createNetCDFFileForModelOutput([runDIR,'anomaly.nc'],shouldOverwriteExisting=1);
model.addNetCDFOutputVariables('waveEnergy')
model.addNetCDFOutputVariables('geostrophicEnergy')

model.integrateToTime(864000/10);

vncfile = model.ncfile;
[A0_i,A0_r,Ap_i,Ap_r,Am_i,Am_r] = vncfile.readVariables('A0_imagp','A0_realp','Ap_imagp','Ap_realp','Am_imagp','Am_realp');
[xFloatT,yFloatT,zFloatT,floatLinPV,floatAPV,time] = vncfile.readVariables('float-x','float-y','float-z','float-LinearPV','float-AvailPV_early','t');
%keyboard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% energy breakup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IW_energy = wvt.waveEnergy;
%VM_energy = wvt.geostrophicEnergy; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot floats, PV flavors, and IC for few floats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xFloatT = xFloatT.';
yFloatT = yFloatT.';
zFloatT = zFloatT.';
floatLinPV = floatLinPV.';
floatAPV = floatAPV.';

figure(1),clf, plot(xFloatT,zFloatT)
xlabel('xFloat')
ylabel('zFloat')

figure(2),clf, plot(xFloatT,yFloatT)
xlabel('xFloat')
ylabel('yFloat')

figure(3),clf, subplot(1,2,1)
plot(time(1:10:end)/inert_per,floatLinPV(1:10:end,56)/(N0^2*f))
xlabel('time (inertial periods)')
ylabel('linPV/N^2/f (unitless)')
title('particle #56')

subplot(1,2,2)
plot(time(1:10:end)/inert_per,floatAPV(1:10:end,56)/f)
xlabel('time (inertial periods)')
ylabel('APV/f (unitless)')
title('particle #56')

figure(5),clf, subplot(1,2,1)
plot(time(1:10:end)/inert_per,floatLinPV(1:10:end,6)/(N0^2)/f)
xlabel('time (inertial periods)')
ylabel('linPV/N^2/f (unitless)')
title('particle #6')

subplot(1,2,2)
plot(time(1:10:end)/inert_per,floatAPV(1:10:end,6)/f)
xlabel('time (inertial periods)')
ylabel('APV/f (unitless)')
title('particle #6')

figure(7),clf, imagesc(x,y,s1(:,:,(Nz-1)/2)')
hold on 
scatter(xFloatT(1,56),yFloatT(1,56),'r','filled')
scatter(xFloatT(1,6),yFloatT(1,6),'k','filled')
colorbar
xlim([18878 24081])
xlabel('x (m)')
ylabel('y (m)')
title('\rho^{prime} (kg/m^3)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set some colors for spectra and limits
colorset = colormap(jet(length(xFloatT)));

xlims = Lx/2+Ly/2*[-1 1];
ylims = [0 Ly];
zlims = Lz/2+Lz/4*[-1 1];
%
for n = 1:length(xFloatT)
  Ap = Ap_r(:,:,:,n)+Ap_i(:,:,:,n);                     % get amplitudes from imag and real parts
  Am = Am_r(:,:,:,n)+Am_i(:,:,:,n);                                     % get amplitudes from imag and real parts
  A0 = A0_r(:,:,:,n)+A0_i(:,:,:,n);                                     % get amplitudes from imag and real parts

  Ap2 = Ap .* conj(Ap);                                         % includes HKE, VKE and PE
  Am2 = Am .* conj(Am);                                 % includes HKE, VKE and PE
  B2 = wvt.A0_TE_factor .* (A0 .* conj(A0));    % no VKE in geostrophic mode

  waveTE_full = 1/Lz * (wvt.Apm_TE_factor .*(Ap2 + Am2)); % full 3D wave total energy
  vortexTE_full = 1/Lz * B2; % full 3D vortex total energy !! THIS DOES INCLUDE BAROTROPIC MODE
  %keyboard
  waveTE_kl = sum(waveTE_full,3); % 2D wave energy (kl)
  vortexTE_kl = sum(vortexTE_full,3); % 2D vortex energy (kl)  !! THIS INCLUDES BAROTROPIC MODE

  waveTE_m(:,n) = squeeze(sum(waveTE_full,[2,1])); % 1D wave energy (m)
  vortexTE_m(:,n) = squeeze(sum(vortexTE_full,[2 1])); % 1D vortex energy (m)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % and set the horizontal wavenumber arrays
  kwavplot = abs(wvt.k(1:Nx/2+1));
  lwavplot = abs(wvt.l(1:Ny/2+1));
  Kh = wvt.Kh;

  % and set the vertical wavenumber array
  wavem = wvt.j/2 * 2*pi/Lz;                                    % (rad/m)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % last step for horizontal wavenumber spectra is to make the horizontal energy arrays into isotropic wavenumber spectra
  % will use k-wavenumber as my isotropic wavenumber
  dk = 2*pi/Lx;
  for nn = 1:length(kwavplot)
    % get indices to wavenumbers within this dk
    ind = find( abs(Kh(:,:,1)) >= (kwavplot(nn)-dk/2) ...
      & abs(Kh(:,:,1)) < (kwavplot(nn)+dk/2) );
    EA_isotr(nn,n) = sum(waveTE_kl(ind));
    EG_isotr(nn,n) = sum(vortexTE_kl(ind));
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure(40);
  % don't clear - need successive spectra to over plot
  if n==1
    clf
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,1)

  temp=2*pi./kwavplot.*EA_isotr(:,n);
  temp(temp==0)=1e-100;
  temp(isnan(temp))=1e-100;
  if n==1
    temp2a = temp;
  end
  loglog(kwavplot/(2*pi),temp,'Color',colorset(n,:),'LineWidth',1)
  ylabel('Horizontal [m^2 s^{-2}] / [cycle/m]')
  hold on
  loglog(kwavplot/(2*pi),temp2a,'Color','b','LineWidth',1,'LineStyle','--')
  hold on
  
 
  colormap(jet)
  title('Wave Energy')
  xlim([1e-5 1e-2])
  ylim([1e-50 1e2])
  grid on
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,2)

  temp=2*pi./kwavplot.*EG_isotr(:,n);
  temp(temp==0)=1e-100;
  temp(isnan(temp))=1e-100;
  if n==1
    temp2b = temp;
  end
  loglog(kwavplot/(2*pi),temp,'Color',colorset(n,:),'LineWidth',1)
  hold on
  loglog(kwavplot/(2*pi),temp2b,'Color','b','LineWidth',1,'LineStyle','--')
  hold on
  
  
  colormap(jet)
  title('Vortex Energy')
  xlim([1e-5 1e-2])
  ylim([1e-50 1e2])
  grid on
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,3)

  temp=2*pi./wavem.*(waveTE_m(:,n));
  temp(temp==0)=1e-100;
  temp(isnan(temp))=1e-100;
  if n==1
    temp2c = temp;
  end
  loglog(wavem/(2*pi),temp,'Color',colorset(n,:),'LineWidth',1)
  ylabel('Vertical [m^2 s^{-2}]/[cycle/m]')
  hold on
  loglog(wavem/(2*pi),temp2c,'Color','b','LineWidth',1,'LineStyle','--')

  colormap(jet)
  xlabel('Wavenumber [cycles/m]')
  xlim([0.5e-3 1e-0])
  ylim([1e-50 1e2])
  grid on
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subplot(2,2,4)

  temp=2*pi./wavem.*(vortexTE_m(:,n));
  temp(temp==0)=1e-100;
  temp(isnan(temp))=1e-100;
  if n==1
    temp2d = temp;
  end
  loglog(wavem/(2*pi),temp,'Color',colorset(n,:),'LineWidth',1)
  hold on
  loglog(wavem/(2*pi),temp2d,'Color','b','LineWidth',1,'LineStyle','--')

  colormap(jet)
  xlabel('Wavenumber [cycles/m]')
  xlim([0.5e-3 1e-0])
  ylim([1e-50 1e2])

  grid on

end