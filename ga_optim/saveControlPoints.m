function saveControlPoints(controlPoints, filename)
    
    fmt = '%10.6f';
    fileID = fopen(filename, 'w');
    for i = 1:numel(controlPoints)
        fprintf(fileID, fmt, controlPoints(i));
        fprintf(fileID, '\n');
    end
    fclose(fileID);
end