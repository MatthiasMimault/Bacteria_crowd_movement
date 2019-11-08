function postLogGeneration(path, caseName, Domain, Space, Attractant, Source, InitBact)
global T R nx CFL typeAt typeInit typeSrc typeObs typeVel A C epsilon...
    BactValue 
% Name = 'TestTest-Test00';
% T = 1000; A = 0.03; C = 0.01; R = 2; eps = 0.5; BactValue = 10;
% nx = 800; CFL = 0.9; 
% typeAt = 'un'; typeSrc = 'un'; typeObs = 'un'; typeVel = 'un';
% typeInit = 'un';
% Domain = [0 0 0 0]; Space = [0 0 0 0]; Attractant = [0 0 0 0]; 
% Source = [0 0 0 0]; InitBact = [0 0 0 0];

fileID = fopen([path '-' caseName '-logs.txt'],'w');
fprintf(fileID, 'Parameters\n');
fprintf(fileID, ['ID Test: ' caseName '\n']);
fprintf(fileID, 'T: %d (s, Maximal time)\n', T);
fprintf(fileID, 'A: %.3f (mm.s-1, Average bacteria velocity)\n', A);
fprintf(fileID, 'C: %.3f (mm2.s-1, Diffusion coefficient)\n', C);
fprintf(fileID, 'R: %.1f (mm, Interaction radius)\n', R);
fprintf(fileID, 'epsilon: %.1f (NU, Interaction strength)\n', epsilon);
fprintf(fileID, 'BactValue: %.1f (cell.mm-2, Cell density)\n', BactValue);
fprintf(fileID, 'nx: %.1f (NU, cell number in X direction)\n', nx);
fprintf(fileID, 'CFL: %.1f (NU, Courant Friedrichs Lewy condition)\n', CFL);
fprintf(fileID, 'typeAt: %s (Attractant type)\n', typeAt);
fprintf(fileID, 'typeSrc: %s (Source type)\n', typeSrc);
fprintf(fileID, 'typeObs: %s (Obstacle type)\n', typeObs);
fprintf(fileID, 'typeVel: %s (Velocity type)\n', typeVel);
fprintf(fileID, 'typeInit: %s (Bacteria initialisation type)\n', typeInit);
fprintf(fileID, 'Domain: [%.2f %.2f %.2f %.2f] (mm, Computational domain verteces)\n'...
    , Domain(1), Domain(2), Domain(3), Domain(4));
fprintf(fileID, 'Space: [%.2f %.2f %.2f %.2f] (mm, Bacteria domain verteces)\n'...
    , Space(1), Space(2), Space(3), Space(4));
fprintf(fileID, 'Attractant: [%.2f %.2f %.2f %.2f] (mm, Attractant region verteces)\n'...
    , Attractant(1), Attractant(2), Attractant(3), Attractant(4));
fprintf(fileID, 'Source: [%.2f %.2f %.2f %.2f] (mm, Source region verteces)\n'...
    , Source(1), Source(2), Source(3), Source(4));
fprintf(fileID, 'InitBact: [%.2f %.2f %.2f %.2f] (mm, Bacteria initialisation region verteces)\n'...
    , InitBact(1), InitBact(2), InitBact(3), InitBact(4));
fclose(fileID);
end

