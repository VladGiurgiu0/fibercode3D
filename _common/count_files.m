function FileCount = count_files(folderPath,ending)
    % Check if the folder exists
    if ~isfolder(folderPath)
        FileCount=NaN;
    end
    
    % Get a list of all files in the folder with .mat extension
    Files = dir(fullfile(folderPath, strcat('*',ending)));
    
    % Count the number of .mat files
    FileCount = length(Files);
end