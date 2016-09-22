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
InputData.Wavelength        = 250:1:900; %nm
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

% InputData.Materials         = {1; 'Au_Jon'; 'glass(BK7)'};
% InputData.Materials         = {1; 'Ag_Jon'; 'glass(BK7)'};
% InputData.Materials         = {1; 'Au_Jon'; 'Ag_Jon'; 'glass(BK7)'};
InputData.Materials         = {1; 'Pd_Plk'; 'YH3'; 'Si_JAW'};

% angle of incidence
% InputData.Theta1 = 10:5:75; % degrees
% InputData.Theta1 = 0; % degrees



 for n_mat=length(InputData.Materials):-1:1
    n(:,n_mat) = interpolate_nk(InputData.Wavelength, InputData.Materials{n_mat}, 1);
    
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
    [ R_p(nwave), T_p(nwave), r_p(nwave)] = films_TMM(InputData.Wavelength(nwave), Epsilon(nwave,:), [0 1]);% last array - polarization fraction: [te tm] 
    [ R_s(nwave), T_s(nwave), r_s(nwave)] = films_TMM(InputData.Wavelength(nwave), Epsilon(nwave,:), [1 0]);
 end
 
 figure
 plot(InputData.Wavelength, R_p, 'Color', 'black'), hold on
 plot(InputData.Wavelength, T_p,  '--', 'Color', 'black'), hold on
 plot(InputData.Wavelength, 1 - R_p - T_p, 'Color', 'red')
 xlabel('Wavelength, nm'), ylabel('R, T')
 legend('R','T','1-R-T')
axis tight

figure
plot(InputData.Wavelength, real(r_p./r_s), '-k', InputData.Wavelength, imag(r_p./r_s), '--k', 'LineWidth', 1.5)
 legend('Re(\rho)','Im(\rho)')
 title(' 4nm Pd / 70 nm Y / Si')
  xlabel('Wavelength, nm'), ylabel('r_p / r_s')
% figure
% plot(InputData.Wavelength, real(r), 'Color', 'black'),hold on
% plot(InputData.Wavelength, imag(r), '--', 'Color', 'black'),
%  xlabel('Wavelength, nm'), ylabel('Re_r, Im_r')
%  legend('Re','Im')
toc;

return
%}

%% Reflection and transmittion dependance on polar angle
%{
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
%}














