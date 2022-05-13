function va_results = voice_analysis_directory(directory)

% Note: We can change here to instead recursively find files, if we want.
fullname = fullfile(directory, '/*.wav');
files = dir(fullname);

result_map = {};

for f = 1:length(files)
    file_name = files(f).name;
    complete_path = fullfile(files(f).folder, file_name);
    [measures_vector, measures_names, F0] = voice_analysis(complete_path);
    result_map{f} = {complete_path, measures_vector, measures_names, F0};    
end

va_results = result_map;

if isdeployed
    save voice_analysis_directory_output.mat va_results
end

end