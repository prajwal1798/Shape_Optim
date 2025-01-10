function exists = checkFilesExist(filePaths)
    
    exists = all(cellfun(@(path) isfile(path), filePaths));
    if exists
        fprintf('All required files are present.\n');
    else
        fprintf('Error: One or more required files are missing.\n');
    end
end