%% Inputs for the bacteria model
global T R nx CFL typeAt typeInit typeSrc typeObs typeVel A D kappa...
    BactValue dev Ro ro

% Main parameters
% T is time in hours, A is bacteria average velocity in mm.h-1, D is diffusion 
% coeff in mm.h-2, R is interaction radius in mm, kappa is interaction 
% strength (0-1), BactValue is in millions of bacteria per mL
T = 5; A = 0.2; D = 0.01; 
R = 0.5; kappa = 0.05;
BactValue = 10;

% Numerics
nx = 400;
Nfiles = 10; 
CFL = 1;
dev = 'dev'; % stale, dev, release

% Domain
% particles - porosity 65%
Ro = 1.18; ro = 1;
Domain = [-8, 8, -10.5, 6];
Space = [-6, 6, -9.2, 5];
Attractant = [-0.2 0.2 4 6]; %[-20 20 20 25] [-1 1 10 25]
InitBact = [-6, 6, -9.2, 3.5];
Source = [-7, 7 -10.2 0];

typeAt = 'root'; %up, root
typeSrc = 'Ulow'; %none, bottom, Ulow
typeObs = 'none'; % none, particles
typeVel = 'att'; % static, att, att2, att-src, none, adv-src
typeInit = 'square'; % none square