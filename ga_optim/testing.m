% clear all 
% close all
% clc

uwallData = load('C:\\optimisation\\uwall.dat');
lwallData = load('C:\\optimisation\\lwall.dat');

[Cl, Cd] = computeLiftDrag(uwallData, lwallData);

fprintf('Lift Coefficient (Cl): %.6f\n', Cl);
fprintf('Drag Coefficient (Cd): %.6f\n', Cd);


