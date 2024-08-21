% This script shows the difference between randomized and maxmin landmark
% read in a variable from a mat file, apply maxmin selector and save to a
% file

clc; clear; close all;
import edu.stanford.math.plex4.*;
% 
% % %filename = '../pentagonsamplesSmall2.mat';

% 
% inputdir = '/Users/nan/Documents/PhD/research/TDA/gudhi/ellipsoids/meshes_nonrigid3d';
% outputdir = '/Users/nan/Documents/PhD/research/TDA/gudhi/ellipsoids/meshes_nonrigid3d/downsampled';
% 
% S = dir(fullfile(inputdir,'*.mat'));
% % S = dir('/Users/nan/Documents/PhD/research/TDA/gudhi/ellipsoids/meshes_nonrigid3d/*.mat');
% nPts = 300;
% seed = 0;
% rng(seed)
% 
% 
% for k = 1:numel(S)
% 
%     filename = fullfile(inputdir,S(k).name);
%     disp(filename)
%     point_cloud = load(filename);
% 
%     fprintf('Processing %s...\n', filename);
% 
%     maxmin_selector = api.Plex4.createMaxMinSelector(point_cloud, nPts);
% 
%     % extract the subset of landmark points from the original point cloud
%     % Note: we need to increment the indices by 1 since Java uses 0-based
%     % arrays
%     maxmin_points = point_cloud(maxmin_selector.getLandmarkPoints() + 1, :);
% 
%     filenameSaved = fullfile(outputdir, sprintf('%s_nPts=%d_seed=%d.mat',filename,nPts,seed));
%     save(filenameSaved, 'maxmin_points')
%     disp(' Done.')
% end

% initialize the point cloud
% filename = '../data/pointsCycloOctane.mat';
filename = 'cyclooctane.mat';
point_cloud = load(filename);
%point_cloud = point_cloud.pentagonsamples;
% point_cloud = point_cloud.pointsCycloOctane;

addPoints = true;

if addPoints == true
    for seed = 0:2
        for nPts = [1000, 800, 500, 300, 100]

            nPtsMax = max(nPts);
    
            rng(seed)
    
            fprintf('Calculating maxmin for seed=%d and nPts=%d...',seed,nPts)
    
            % create the landmark selectors
            maxmin_selector = api.Plex4.createMaxMinSelector(point_cloud, nPts);
    
            % extract the subset of landmark points from the original point cloud
            % Note: we need to increment the indices by 1 since Java uses 0-based
            % arrays
            maxmin_points = point_cloud(maxmin_selector.getLandmarkPoints() + 1, :);
    
            filenameSaved = sprintf('%s_nMax=%d_nPts=%d_seed=%d.mat',filename,nPtsMax,nPts,seed);
            save(filenameSaved, 'maxmin_points')
            disp(' Done.')

            point_cloud = maxmin_points;
        end
    end

else

    for seed = 0:2
        for nPts = [100, 300, 500, 800, 1000]
    
            rng(seed)
    
            fprintf('Calculating maxmin for seed=%d and nPts=%d...',seed,nPts)
    
            % create the landmark selectors
            maxmin_selector = api.Plex4.createMaxMinSelector(point_cloud, nPts);
    
            % extract the subset of landmark points from the original point cloud
            % Note: we need to increment the indices by 1 since Java uses 0-based
            % arrays
            maxmin_points = point_cloud(maxmin_selector.getLandmarkPoints() + 1, :);
    
            filenameSaved = sprintf('%s_nPts=%d_seed=%d.mat',filename,nPts,seed);
            save(filenameSaved, 'maxmin_points')
            disp(' Done.')
        end
    end

end
