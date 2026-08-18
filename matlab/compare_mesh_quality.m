% MATLAB Script: compare_mesh_quality.m
% Reads HDF5 mesh quality files generated before and after optimization
% and plots comparative histograms for key hexahedral quality metrics.

clear; clc; close all;

% Search for generated HDF5 quality files
files_before = dir('mesh_before_opt_*.h5');
files_after  = dir('mesh_after_opt_*.h5');

if isempty(files_before)
    files_before = dir('../mesh_before_opt_*.h5');
    files_after  = dir('../mesh_after_opt_*.h5');
end

if isempty(files_before) || isempty(files_after)
    error('Could not find mesh_before_opt_*.h5 or mesh_after_opt_*.h5 files. Please run hexmesh first.');
end

file_before = fullfile(files_before(1).folder, files_before(1).name);
file_after  = fullfile(files_after(1).folder, files_after(1).name);

fprintf('=========================================================\n');
fprintf(' HexMesh Quality Comparison (MATLAB)\n');
fprintf(' Before: %s\n', file_before);
fprintf(' After:  %s\n', file_after);
fprintf('=========================================================\n\n');

% Metrics to load and compare
metrics = { ...
    'ScaledJacobian',   'Scaled Jacobian [-1, 1] (Ideal: 1.0)'; ...
    'ConditionNumber',  'Condition Number [1, inf) (Ideal: 1.0)'; ...
    'Skew',             'Principal Axis Skew [0, 1] (Ideal: 0.0)'; ...
    'Shape',            'Shape Metric (0, 1] (Ideal: 1.0)'; ...
    'MinFaceAngle',     'Min Face Angle (Degrees) (Ideal: 90.0)'; ...
    'EdgeRatio',        'Edge Ratio [1, inf) (Ideal: 1.0)' ...
};

figure('Name', 'HexMesh Quality Comparison: Before vs After Optimization', 'Position', [100, 100, 1200, 800]);

for i = 1:size(metrics, 1)
    dataset_name = ['/Sem3D/', metrics{i,1}];
    
    try
        data_before = h5read(file_before, dataset_name);
        data_after  = h5read(file_after, dataset_name);
    catch ME
        warning('Could not read dataset %s: %s', dataset_name, ME.message);
        continue;
    end
    
    subplot(2, 3, i);
    
    % Cap extreme condition number / edge ratio for visualization if needed
    if strcmp(metrics{i,1}, 'ConditionNumber') || strcmp(metrics{i,1}, 'EdgeRatio')
        cap_val = 50;
        data_before(data_before > cap_val) = cap_val;
        data_after(data_after > cap_val)   = cap_val;
    end
    
    histogram(data_before, 40, 'FaceColor', [0.2 0.4 0.8], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Before Opt');
    hold on;
    histogram(data_after, 40, 'FaceColor', [0.8 0.3 0.1], 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'After Opt');
    hold off;
    
    title(metrics{i,2}, 'FontSize', 11, 'FontWeight', 'bold');
    xlabel('Metric Value');
    ylabel('Element Count');
    legend('Location', 'best');
    grid on;
    
    % Print statistics summary
    fprintf('Metric: %s\n', metrics{i,1});
    fprintf('  BEFORE: min = %8.4f, mean = %8.4f, 5th percentile = %8.4f\n', ...
        min(data_before), mean(data_before), prctile(data_before, 5));
    fprintf('  AFTER:  min = %8.4f, mean = %8.4f, 5th percentile = %8.4f\n\n', ...
        min(data_after), mean(data_after), prctile(data_after, 5));
end

sgtitle('HexMesh Mesh Quality Distribution Before vs. After Optimization', 'FontSize', 14, 'FontWeight', 'bold');
