function [domAt, domBd, domDef, domSrc] ...
    = fRegionGeneration(X, Domain, Space, Attractant, Source)
% Define the numerical domain of definition, attractant, boundary and
% source based on geometric input and a grid X;

domAt = zeros(size(X));
domBd = zeros(size(X));
domDef = zeros(size(X));
domSrc = zeros(size(X));