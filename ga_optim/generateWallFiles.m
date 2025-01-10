function [uwallData, lwallData] = generateWallFiles(controlPoints)

    % INITIALIZE NURBS FUNCTION 
    fileNameIGES = 'naca0012.igs';
    nsd = 2;
    quadRef = defineQuadratureAdaptive();
    nurbs = nurbsReadIgesBoundary(fileNameIGES, nsd);
    nurbs = nurbsSetupStruct(nurbs, quadRef, nsd);

    % LOAD UPDATED CONTROLPOINTS
    controlPoints = load('C:\\optimisation\\controlPoints.dat');

    % UPDATE NURBS(4) 
    nurbs(4).Pw(2:7, 1) = controlPoints(1:2:11);  
    nurbs(4).Pw(2:7, 2) = controlPoints(2:2:12);

    % UPDATE NURBS(5)   
    nurbs(5).Pw(2:6, 1) = controlPoints(13:2:21); 
    nurbs(5).Pw(2:6, 2) = controlPoints(14:2:22); 

    nu = 0.5;
    nurbs(5).Pw(7, 1) = nurbs(5).Pw(8, 1) + nu * (nurbs(5).Pw(8, 1) - nurbs(4).Pw(7, 1));
    nurbs(5).Pw(7, 2) = nurbs(5).Pw(8, 2) + nu * (nurbs(5).Pw(8, 2) - nurbs(4).Pw(7, 2));

    nurbs = nurbsSetupStruct(nurbs, quadRef, nsd);

    % INITIALIZE WALL DATA
    uwallData = zeros(140, 7);
    lwallData = zeros(140, 7);

    % COMPUTE GEO DATA [COORDS/NORMALS]   
    nPoints = 140;
    uSamples = linspace(0, 5, nPoints);

    % UPPER SURFACE GEO DATA
    for i = 1:nPoints
        [pt, dpt] = nurbsCurveDerivPoint(nurbs(4), uSamples(i));
        normal = [-dpt(2), dpt(1)]; 
        normal = normal / norm(normal);  
        uwallData(i, 1:2) = pt(1:2);
        uwallData(i, 6:7) = normal;
    end

    % LOWER SURFACE GEO DATA
    for i = 1:nPoints
        [pt, dpt] = nurbsCurveDerivPoint(nurbs(5), uSamples(i));
        normal = [dpt(2), -dpt(1)];  
        normal = normal / norm(normal);
        lwallData(i, 1:2) = pt(1:2); 
        lwallData(i, 6:7) = normal;
    end

    % fprintf('%.6f \n' , controlPoints);
end


