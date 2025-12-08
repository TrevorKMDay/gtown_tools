function call_LI(files, roi, thr1, output_file, overwrite, spm)

% This sets default values
arguments
    files
    roi
    thr1
    output_file
    overwrite = 0
    spm = '/Users/tkmd/Documents/MATLAB/spm12'
end

% File input ====
files_split = split(files) ;
files_split(cellfun("isempty", files_split)) = [];
files_split_char = char(files_split) ;

disp(files_split_char)

n_files = size(files_split_char) ;
disp(strcat("[call_LI] I was given ", string(n_files(1)), " files."))

% overwrite = 0     don't do anything if output file exists
%           = 1     delete output file and start over
%           = 2     append new lines to existing file (not tested)

% TO DO: Check if # lines in output matches n_files

if exist(output_file, "file") && overwrite == 0
    % End execution if file exists
    disp(strcat("[call_LI] Output file ", output_file, ...
                " exists and overwrite not set, not overwriting."))
    return
elseif exist(output_file, "file") && overwrite == 1
    disp(strcat("[call_LI] Overwriting output file ", output_file, "!"))
    delete(output_file)
elseif exist(output_file, "file") && overwrite == 2
    disp(strcat("[call_LI] Appending to ", output_file, "!"))
end

% add libraries
disp(strcat("[call_LI] SPM: ", spm))
addpath(spm);
addpath(strcat(spm, '/toolbox/LI'));

disp(strcat("[call_LI] ROI: ", roi))

% Thresholding ====
% thr1    =  1 (use same threshold for all images)
%         =  0 (individual thresholding for all images)
%         = -1 (adaptive thresholding
%         = -2 (rank-based thresholding)
%         = -3 (iterative thresholding, LI-curves)
%         = -4 (no threshold)
%         = -5 (bootstrapping)

disp(strcat("[call_LI] Using threshold setting: ", thr1))

% B1 = the ROI
% C1 = exclusive mask [1: midline 5mm]

li_input = struct("A", files_split_char, ...
                    "B1", roi, ...
                    "C1", 1, ...
                    "thr1", str2double(thr1), ...
                    "pre", 0, ...
                    "outfile", output_file) ;

LI(li_input)

disp(strcat("[call_LI] Done with ", output_file))

end
