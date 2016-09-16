function [ R, T, r ] = films_TMM( lam0, epsilon, theta, mu  )
% This MATLAB program implements the transfer matrix method
% INITIALIZE MATLAB

% close all;
% clc
% clear all;

tic

InputData.MATLABDir = 'D:\works\';

% commom functions
addpath([InputData.MATLABDir,'_Data'],[InputData.MATLABDir,'_Common'],[InputData.MATLABDir,'_Algorithms'])
% common globals
load commonSI
% UNITS
degrees = pi/180;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DEFINE SIMULATION PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOURCE PARAMETERS
if ~exist('lam0', 'var')
    fprintf('\n Variable lam0 does not exist\n')
%     lam0 = 1;
    lam0 = 633; % nm;  free space wavelength
end	
if ~exist('theta', 'var')
    fprintf('\n Variable THETA does not exist\n')
    theta = 70 * degrees; %elevation angle
else
    theta = theta * degrees;
end
phi = 0 * degrees; %azimuthal angle
pte = 0; %amplitude of TE polarization
ptm = 1; %amplitude of TM polarization
% EXTERNAL MATERIALS
if ~exist('epsilon','var')
    fprintf('\n Variable epsilon does not exist\n')
    er1 = 1.0; %permittivity in the reflection region
    er2 = 2.2954 - 1i*4.8254e-23; %permittivity in the transmission region
%     ER = [ (1.37539 - 1i*1.26657)^2 ];% = 0.2875 - 3.4841i
    ER = [ -11.8183 - 1i*1.2367];
%     ER = [ (1.55)^2 ]; %array of permittivities in each layer
else
    er1 = epsilon(1); %permittivity in the reflection region
    er2 = epsilon(end); %permittivity in the transmission region
    ER  = epsilon(2:end-1); %array of permittivities in each layer
end
if ~exist('mu','var')
%     fprintf('\n Variable MU does not exist\n')
    ur1 = 1.0; %permeability in the reflection region
    ur2 = 1.0; %permeability in the transmission region
    UR = ones(size(ER)); %array of permeability in each layer
else
    ur1 = mu(1); %permeability in the reflection region
    ur2 = mu(end); %permeability in the transmission region
    UR  = mu(2:end-1); %array of permeability in each layer
end
% % LAYERS
% UR = [ 1 ]; %array of permeabilities in each layer
% ER = [ 0.916*(0.916 - 1i*1.84) ]; %array of permittivities in each layer
L = [ 55, 600000 ]; % nm; array of the thicknesses of each layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


kx = sqrt(ur1*er1)*sin(theta)*cos(phi);
ky = sqrt(ur1*er1)*sin(theta)*sin(phi);

k0 = 2*pi/lam0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gap medium parameters   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kz0 = sqrt(1 - kx^2 - ky^2);
if kz0 == 0
    fprintf('\n WARNING! kz0 = 0! \n ')
end
Q0 = [kx*ky, 1 - kx^2; ky^2 - 1, -kx*ky];
OMEGA0 = 1i*kz0*eye([2 2]);
V0 = Q0 / OMEGA0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Initialization of global S-matrix   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S11g = zeros([2 2]);
S12g = eye([2 2]);
S21g = eye([2 2]);
S22g = zeros([2 2]);
t1 = toc;
NLAY = length(L);
for nlay = 1:NLAY
    kzi = sqrt(UR(nlay)*ER(nlay) - kx^2 - ky^2);
    if kzi == 0
        fprintf('\nWARNING! Total reflection in layer %d \n', nlay);
    end
    Qi = [kx*ky, UR(nlay)*ER(nlay) - kx^2; ky^2 - UR(nlay)*ER(nlay), -kx*ky]/UR(nlay);
    OMEGAi = 1i*kzi*eye([2 2]);
    Vi = Qi/OMEGAi;
    
    Ai = eye([2 2]) + Vi\V0;
    Bi = eye([2 2]) - Vi\V0;
    Xi = expm(-OMEGAi*k0*L(nlay));
    
%     S11i = (Ai*Ai - Xi*Bi*Xi*Bi)\(Xi*Bi*Xi*Ai - Ai*Bi); %test
%     S12i = (Ai*Ai - Xi*Bi*Xi*Bi)\(Ai*Xi*Ai - Xi*Bi*Bi); %test
    
    S11i = (Ai - Xi*Bi/Ai*Xi*Bi)\(Xi*Bi/Ai*Xi*Ai - Bi);
    S12i = (Ai - Xi*Bi/Ai*Xi*Bi)\Xi*(Ai - Bi/Ai*Bi);
    S21i = S12i;
    S22i = S11i;
    
%     [c11i,c12i,c21i,c22i] = CHECK_SH_invS(S11i,S12i,S21i,S22i)

    CoefF = S12g/(eye([2 2]) - S11i*S22g);
    CoefD = S21i/(eye([2 2]) - S22g*S11i);
    
    S11g = S11g + CoefF*S11i*S21g;
    S12g = CoefF*S12i;
    S21g = CoefD*S21g;
    S22g = S22i + CoefD*S22g*S12i;
%     [c11g,c12g,c21g,c22g] = CHECK_SH_invS(S11g,S12g,S21g,S22g)

end

t2 = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  S-matrix of reflection region   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kzref = sqrt(ur1*er1 - kx^2 - ky^2);
Qref = [kx*ky, ur1*er1 - kx^2; ky^2 - ur1*er1, -kx*ky]/ur1;
OMEGAref = 1i*kzref*eye([2 2]);
Vref = Qref/OMEGAref;

Aref = eye([2 2]) + V0\Vref;
Bref = eye([2 2]) - V0\Vref;

S11ref = -Aref\Bref;
S12ref = 2*inv(Aref);
S21ref = 0.5*(Aref - Bref/Aref*Bref);
S22ref = Bref/Aref;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S-matrix of transmission region  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kztrn = sqrt(ur2*er2 - kx^2 - ky^2);
Qtrn = [kx*ky, ur2*er2 - kx^2; ky^2 - ur2*er2, -kx*ky]/ur2;
OMEGAtrn = 1i*kztrn*eye([2 2]);
Vtrn = Qtrn/OMEGAtrn;

Atrn = eye([2 2]) + V0\Vtrn;
Btrn = eye([2 2]) - V0\Vtrn;

S11trn = Btrn/Atrn;
S12trn = 0.5*(Atrn - Btrn/Atrn*Btrn);
S21trn = 2*inv(Atrn);
S22trn = -Atrn\Btrn;
% [E11r,E12r,E21r,E22r] = CHECK_SH_S(S11ref,S12ref,S21ref,S22ref)
% [C11r,C12r,C21r,C22r] = CHECK_SH_invS(S11ref,S12ref,S21ref,S22ref)
% [E11t,E12t,E21t,E22t] = CHECK_SH_S(S11trn,S12trn,S21trn,S22trn)
% [C11t,C12t,C21t,C22t] = CHECK_SH_invS(S11trn,S12trn,S21trn,S22trn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Connecting S-matrix of device to outer medium matrixes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CoefD = S12ref/(eye([2 2]) - S11g*S22ref);
CoefF = S21g/(eye([2 2]) - S22ref*S11g);

S22g = S22g + CoefF*S22ref*S12g;
S21g = CoefF*S21ref;
S12g = CoefD*S12g;
S11g = S11ref + CoefD*S11g*S21ref;

% [E11g,E12g,E21g,E22g] = CHECK_SH_S(S11g,S12g,S21g,S22g)
% [C11g,C12g,C21g,C22g] = CHECK_SH_invS(S11g,S12g,S21g,S22g)

CoefM = S12g/(eye([2 2]) - S11trn*S22g);
CoefO = S21trn/(eye([2 2]) - S22g*S11trn);

S11g = S11g + CoefM*S11trn*S21g;
S12g = CoefM*S12trn;
S21g = CoefO*S21g;
S22g = S22trn + CoefO*S22g*S12trn;
% [E11g,E12g,E21g,E22g] = CHECK_SH_S(S11g,S12g,S21g,S22g)
% [C11g,C12g,C21g,C22g] = CHECK_SH_invS(S11g,S12g,S21g,S22g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Source Polarization vector     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if theta == 0
    cinc = [ptm; pte];
else
       P = [ptm*cos(theta)*cos(phi) - pte*sin(phi), ptm*cos(theta)*sin(phi) + pte*cos(phi), -ptm*sin(theta)];
    cinc = [ptm*cos(theta)*cos(phi) - pte*sin(phi); ptm*cos(theta)*sin(phi) + pte*cos(phi)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Transmitted and reflected fields  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Eref = S11g*cinc;
Etrn = S21g*cinc;
%%%%%% Ez:
Ezref = -(kx*Eref(1) + ky*Eref(2))/kzref;
Eztrn = -(kx*Etrn(1) + ky*Etrn(2))/kztrn;

r = Eref(1) + Eref(2) + Ezref

R = abs(Eref(1))^2 + abs(Eref(2))^2 + abs(Ezref)^2
T = (abs(Etrn(1))^2 + abs(Etrn(2))^2 + abs(Eztrn)^2)*real(ur1*kztrn/ur2/kzref)
% fprintf('\n loop time %d \n', t2-t1)