function logControlPoints(controlPoints)
    
    controlPoints = controlPoints(:);

    % UPDATE CONTROLPOINTS DATA 
    save('C:\\optimisation\\controlPoints.dat', 'controlPoints', '-ascii');
end
