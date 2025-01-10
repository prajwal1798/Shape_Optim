function cd_cl_ratio = objectiveFunction(controlPoints)

    % LOG CONTROLPOINTS DATA
    logControlPoints(controlPoints);
    
    % GENERATE GEO INFO FOR WALL DATA
    [uwallData, lwallData] = generateWallFiles(controlPoints);
    
    % SAVE GEO INFO TO WALL DATA
    saveWallData(uwallData, 'C:\\optimisation\\uwall.dat');
    saveWallData(lwallData, 'C:\\optimisation\\lwall.dat');

    % PREDICT AND SAVE WALL DATA
    uwallData = runPrediction(uwallData, 'upper');
    lwallData = runPrediction(lwallData, 'lower');

    % COMPUTE SURFACE INTEGRALS 
    [Cl, Cd] = computeLiftDrag(uwallData, lwallData);
    
    fprintf('Drag Coefficient: %10.6f \n', Cd);
    fprintf('Lift Coefficient: %10.6f \n', Cl);
    
    % COMPUTE OBJ. FUNCTION
    % Check for negative or zero values
    if Cl <= 0 || Cd <= 0
        % Apply a large penalty if lift or drag is non-positive
        cd_cl_ratio = 1e6;  % High penalty value
    else
        % Calculate the drag-to-lift ratio
        cd_cl_ratio = Cd / Cl;
    end
end
