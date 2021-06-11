# Bacteria Crowd Movement model
Copyright (C) 2021  Matthias Mimault
Contact matthias.mimault -at- hutton.ac.uk

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

# Description
Matlab code. 2D finite volume model with nonlocal interaction

Features: Zero flux Neumann boundary, source boundary, attraction (eikonal 
equation), attraction (neighbourhood), diffusion, vector field timer switch, 
log generator

# Usage
To run
Go to Exec and run the desired run script (for instance run.m)
You can modify its parameters in section 0 (Domain, density, interaction 
strenght, diffusion, time and discretisation)

To plot
Go to Exec and run the desired plot script (for instance plot.m)
You can choose a type of output (fig or png) and a scale of density (min-max)
