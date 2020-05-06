%% Script to fix the vector field generation
close all
clear

%% 0. Settings
% Names
caseTitle = 'DBG-VectorField';
caseDate = '0506';
runnb = '1';
suffix = '';
global Dx Dy nx

% input file
input_V34

%% Nothing to be modified below
typeRepulsion = 0;

% Directories
addpath('..\Include')
addpath('..\Source')
addpath('..\Cases')
[dataRoot, caseName, fileName] = fFolderMaker( ...
    caseTitle, caseDate, runnb, nx, suffix);

%% 1. Initialisation grid and dependant parameters
[X,Y,Dx,Dy] = fGridGeneration_dbg(nx,Domain);

[domAt, domBd, domDef, domSrc] ...
    = fRegionGeneration(X, Y, Space, Attractant, Source);

% Generation velocity field
[Vxo, Vyo, PhiC] = fVelocityGeneration_dbg(X, Y, domAt, domBd, domSrc);
Vx = Vxo;
Vy = Vyo;