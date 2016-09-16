% MATLAB my program folder
clear all
InputData.MATLABDir = 'D:\works\';

% commom functions
addpath([InputData.MATLABDir,'_Data'],[InputData.MATLABDir,'_Common'],[InputData.MATLABDir,'_Algorithms'])
% common globals
load commonSI

InputData.OutputDir = '.\';

  
tic
%% Definitions
% Spectra definition
InputData.Wavelength        = 200:1:900; %nm
% InputData.Wavelength        = 633; %nm
% InputData.Wavelength       = 515; %nm

InputData.Polarization_inc = [0 1]; % fraction of [TE TM] polarization

% material parameters
% InputData.Thickness           = 5; %nm, thickness of the metal slab
% dielectric
% InputData.Materials         = {1; 1.37539 - 1i*1.26657; 1};%permittivity
% InputData.Materials         = {1; 0.1796 - 1i*3.4425; 1.5151}; %gold and BK7 at 633 nm
% InputData.Materials         = {1; struct('type', 'maxwell_garnett_eps', 'nk_medium', 'glass(BK7)', 'nk_incl',  'Au_Jon', 'vol_fr', 0); 1};%permittivity
% InputData.Materials         = {1; struct('type', 'maxwell_garnett_eps', 'nk_medium', 'Au_Jon', 'nk_incl',  'Ag_Jon', 'vol_fr', 0.7); 1};%permittivity
% InputData.Materials         = {1; 'Au_Jon'; 1}; %{environment, metalic film, substrate}
% InputData.Materials         = {1; 'Au_Jon'; 'glass(BK7)'}; %{environment, metalic film, substrate}

% InputData.Materials         = {1; 'Au_Jon'; 1};
% InputData.Materials         = {1; 'Ag_Jon'; 1};
% InputData.Materials         = {1; 'glass(BK7)'; 1};
InputData.Materials         = {1; 'Au_Jon'; 'Si_JAW'; 1};

% angle of incidence
% InputData.Theta1 = 10:5:75; % degrees
% InputData.Theta1 = 0; % degrees

t0 = toc;

 for n_mat=length(InputData.Materials):-1:1
    n(:,n_mat) = interpolate_nk(InputData.Wavelength, InputData.Materials{n_mat});
    
 end % (n_0, n_1, n_2)
 Epsilon = n.^2;

 figure
 plot(InputData.Wavelength, real(Epsilon(:,2)), '-k', InputData.Wavelength, imag(Epsilon(:,2)), '-r')
 xlabel('Wavelength, nm'), ylabel('epsilon')
 legend('Re','Im')
 
%% Reflection and transmittion dependance on wavelength
% %{
 N_WAVE=length(InputData.Wavelength);
 for nwave=N_WAVE:-1:1
    [ R(nwave), T(nwave), r(nwave)] = films_TMM(InputData.Wavelength(nwave), Epsilon(nwave,:), InputData);
 end
 
 figure
 plot(InputData.Wavelength, R, 'Color', 'black'), hold on
 plot(InputData.Wavelength, T,  '--', 'Color', 'black'), hold on
 plot(InputData.Wavelength, R+T, 'Color', 'red')
 xlabel('Wavelength, nm'), ylabel('R, T')
 legend('R','T','R+T')
axis tight
% figure
% plot(InputData.Wavelength, real(r), 'Color', 'black'),hold on
% plot(InputData.Wavelength, imag(r), '--', 'Color', 'black'),
%  xlabel('Wavelength, nm'), ylabel('Re_r, Im_r')
%  legend('Re','Im')
 tend = toc;
fprintf('\n total time - %d\n', tend - t0)
return
%}

%% Reflection and transmittion dependance on polar angle
THETA = 0:1:89; %degrees
N_THETA=length(THETA);
for ntheta=N_THETA:-1:1
    [ R(ntheta), T(ntheta), r(ntheta)] = films_TMM(InputData.Wavelength, Epsilon ,THETA(ntheta), Polarization_inc );
end

figure
plot(THETA, R, 'Color', 'black'), hold on
plot(THETA, T, '--', 'Color', 'black')
plot(THETA, 1-R-T, 'Color', 'red')
 xlabel('Theta, degrees'), ylabel('R, T')
 legend('R','T','Absorbance')















