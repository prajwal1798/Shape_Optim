function controlPoints = initializeControlPoints()
    % INITIALIZE NURBS -> NACA0012
    fileNameIGES = 'naca0012.igs';
    nsd = 2;
    nurbs = nurbsReadIgesBoundary(fileNameIGES, nsd);

    % SET BOUNDARY POINTS
    nurbs(4).Pw(1, 1:2) = [0.5, 0.0];    
    nurbs(5).Pw(8, 1:2) = [-0.5, 0.0];   
    nurbs(4).Pw(8, 1:2) = [-0.5, 0.0];   
    nurbs(5).Pw(1, 1:2) = [0.5, 0.0];    
    
    % INITIALIZE CONTROLPOINTS DATA
    controlPoints = zeros(22, 1);

    % POPULATE WITH NACA0012 COORDINATES (UPPER SURFACE)
    controlPoints(1:2:11) = nurbs(4).Pw(2:7, 1);  
    controlPoints(2:2:12) = nurbs(4).Pw(2:7, 2);  

    % POPULATE WITH NACA0012 COORDINATES (LOWER SURFACE)
    controlPoints(13:2:21) = nurbs(5).Pw(2:6, 1); 
    controlPoints(14:2:22) = nurbs(5).Pw(2:6, 2); 
    
    logControlPoints(controlPoints);
end
