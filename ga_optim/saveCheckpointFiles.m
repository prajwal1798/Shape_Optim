function saveCheckpointFiles(uwallData, lwallData, iteration)
    
    upperFile = sprintf('C:\\optimisation\\Checkpoints\\uwall_checkpoint_%d.dat', iteration);
    lowerFile = sprintf('C:\\optimisation\\Checkpoints\\lwall_checkpoint_%d.dat', iteration);
    
    saveWallData(uwallData, upperFile);
    saveWallData(lwallData, lowerFile);
end

function saveOptimalSolution(controlPoints, uwallData, lwallData)
    
    saveWallData(uwallData, 'C:\\optimisation\\uwall_optimal.dat');
    saveWallData(lwallData, 'C:\\optimisation\\lwall_optimal.dat');
    saveControlPoints(controlPoints, 'C:\\optimisation\\optimal_controlpoints.dat');
end

function saveWallData(data, filename)
    
    fmt = '%18.10e';
    fileID = fopen(filename, 'w');
    for i = 1:size(data, 1)
        fprintf(fileID, fmt, data(i, :));
        fprintf(fileID, '\n');
    end
    fclose(fileID);
end

function saveControlPoints(controlPoints, filename)
    
    save(filename, 'controlPoints', '-ascii');
end
