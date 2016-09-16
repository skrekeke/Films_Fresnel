% MATLAB my program folder
clear all
InputData.MATLABDir = 'D:\works\';

% commom functions
addpath([InputData.MATLABDir,'_Data'],[InputData.MATLABDir,'_Common'],[InputData.MATLABDir,'_Algorithms'])
% common globals
load commonSI

InputData.OutputDir = '.\';



%% Definitions
% Spectra definition
InputData.Wavelength        = 400; %nm
% InputData.Wavelength       = 515; %nm

% material parameters
% InputData.Thickness           = 5; %nm, thickness of the metal slab
% dielectric
InputData.Materials         = {1; 3; 1.5};%test R+T=1
% InputData.Materials         = {1; 'Au_Jon'; 1};
% InputData.Materials         = {1; 'Au_Jon'; 1.55}; %{environment, metalic film, substrate}

% angle of incidence
% InputData.Theta1 = 10:5:75; % degrees
InputData.Theta1 = 0; % degrees
%% Reflectance and transmittance Coeficients of metal slab on dielectric substrate
%(s, p polarizations)

% %{
 for n_mat=length(InputData.Materials):-1:1
    n(:,n_mat) = interpolate_nk(InputData.Wavelength, InputData.Materials{n_mat});
    N(:,n_mat) = sqrt(real(n(:,n_mat)).*(real(n(:,n_mat)) - 1i*imag(n(:,n_mat))));
 end % (n_0, n_1, n_2)

% [R_s, T_s,~,~] = RT_metalic_slab(InputData)

% [R_s, T_s] = RT_metalic_slab_2(InputData)

Theta = [repmat(InputData.Theta1, length(InputData.Wavelength), 1), asind(n(:,1).*sind(InputData.Theta1)./n(:,2)),asind(n(:,1).*sind(InputData.Theta1)./n(:,3))] %(Theta0, Theta1, Theta2)
% alpha =  2*pi*InputData.Thickness*cosd(Theta(:,2))./InputData.Wavelength' % delt0, delta1 thickness [nm]/ wavelength [nm]
% 
% 
% for n_spectra=length(InputData.Wavelength):-1:1
%     [R_s(n_spectra), T_s(n_spectra), DMs(n_spectra),DMp(n_spectra)] = RT_metalic_slab_M(alpha(n_spectra), Theta(n_spectra,:),n(n_spectra,:));
% end
Thicknesses = [5:5:450]
for thick=length(Thicknesses):-1:1
    InputData.Thickness = Thicknesses(thick)
    [R_s(thick), T_s(thick), R_p(thick), T_p(thick)] = RT_metalic_slab(InputData)

%     [R_s(thick), T_s(thick)] = RT_metalic_slab_2(InputData)
    
%     alph = 2*pi*InputData.Thickness*cosd(Theta(:,2))./InputData.Wavelength'
%     [R_s(thick), T_s(thick), DMs(thick),DMp(thick)] = RT_metalic_slab_M(alph, Theta, N)
end    

figure
plot(Thicknesses.*N(2)./InputData.Wavelength, R_s), hold on
plot(Thicknesses.*N(2)./InputData.Wavelength, T_s),
plot(Thicknesses.*N(2)./InputData.Wavelength, T_s + R_s),
xlabel('Thickness*n, lambda'), ylabel('R, T (TE)'),
legend('R_s','T_s', 'R_s + T_s')

figure,
% plot(InputData.Wavelength, real(r_p)), xlabel('Wavelength, nm'), ylabel('Reflection coefficient'), hold on %
plot(InputData.Wavelength, R_s), hold on
plot(InputData.Wavelength, T_s),
plot(InputData.Wavelength, T_s + R_s),
xlabel('Wavelength, nm'), ylabel('R, T (TE)'),
legend('R_s','T_s', 'R_s + T_s')

figure,
plot(InputData.Wavelength, R_p), hold on
plot(InputData.Wavelength, T_p),
plot(InputData.Wavelength, T_p + R_p),
xlabel('Wavelength, nm'), ylabel('R, T (TM)'),
legend('R_p', 'T_p', 'R_p + T_p')
%}
%% Radial dependence
%{
InputData.Theta = 0:10:90;
for n_angle=length(InputData.Theta):-1:1
    InputData.Theta1 = InputData.Theta(n_angle)
    [R_s(:,n_angle), T_s(:,n_angle)] = RT_metalic_slab(InputData)
end
for wavelength=length(InputData.Wavelength):-1:1
    
    figure,
    plot(InputData.Theta, R_s(wavelength,:)), hold on
    plot(InputData.Theta, T_s(wavelength,:)),
    plot(InputData.Theta, T_s(wavelength,:) + R_s(wavelength,:)),
    xlabel('Theta1, degrees'), ylabel('R, T (TE)'),
    legend('R_s','T_s', 'R_s + T_s','Lambda = %f nm', wavelength)
end
%}