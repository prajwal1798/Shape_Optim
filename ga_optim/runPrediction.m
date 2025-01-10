function predictedValues = runPrediction(wallData, surfaceType)
    pyExec = 'C:\Users\prajw\AppData\Local\Programs\Python\Python312\python.exe';
    pyScript = 'C:\\optimisation\\predictWallValues.py';
    
    % SET TEXT INTERPRETER
    command = sprintf('set PYTHONIOENCODING=UTF-8 && "%s" "%s" "%s"', pyExec, pyScript, surfaceType);
    
    % CALL .py SCRIPT TO PREDICT
    [status, cmdout] = system(command);

    % CHECK
    if status == 0
        disp('Python script executed successfully');
    else
        disp(['Python script execution failed: ', cmdout]);
    end

    % LAOD PREDICTIONS FROM .py & SAVE TO WALL DATA
    if strcmp(surfaceType, 'upper')
        predictedData = load('C:\\optimisation\\uwall.dat');
    else
        predictedData = load('C:\\optimisation\\lwall.dat');
    end
    
    % UPDATE WALL DATA WITH PREDICTIONS
    wallData(:, 3:5) = predictedData(:, 3:5);
    predictedValues = wallData;    
end
