function saveWallData(data, filename)    
    fmt = '%10.6f';
    fileID = fopen(filename, 'w');
    for i = 1:size(data, 1)
        fprintf(fileID, fmt, data(i, :));
        fprintf(fileID, '\n');
    end
    fclose(fileID);
end