%% Inputs for the bacteria model
% Top attraction, source at bottom
global T R nx CFL typeAt typeInit typeSrc typeObs typeVel A D kappa...
    BactValue dev Ro ro rateCol


% Main parameters
% T is time in hours, A is bacteria average velocity in mm.h-1, D is diffusion 
% coeff in mm.h-2, R is interaction radius in mm, kappa is interaction 
% strength (0-1), BactValue is in millions of bacteria per mL
T = 15; A = 0.2; D = 0.025; 
R = 0.1; kappa = 0.1;
BactValue = 15; rateCol = 0.2;

% Numerics
nx = 400;
Nfiles = 20; 
CFL = 1;
dev = 'release'; % stale, dev, release

% Domain
% particles - porosity 65%
Ro = 1.18; ro = 1;
Domain = [-8, 8, -10.5, 6];
Space = [-6, 6, -9.2, 5];
Attractant = [-0.2 0.2 4 6]; %[-20 20 20 25] [-1 1 10 25]
InitBact = [-6, 6, -9.2, 5];
Source = [0, 0, 0, 0];

typeAt = 'root'; %up, root
typeSrc = 'none'; %none, bottom, Ulow
typeObs = 'none'; % none, particles
typeVel = 'att'; % static, att, att2, att-src, none, adv-src
typeInit = 'square'; % none square