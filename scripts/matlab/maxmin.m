% Downsample a dataset so that the points are most "spread-out".
% See https://appliedtopology.github.io/javaplex/doc/edu/stanford/math/plex4/metric/landmark/MaxMinLandmarkSelector.html
% for more info

clc; clear; close all;
import edu.stanford.math.plex4.*;

filename = 'cyclooctane.mat';
point_cloud = load(filename);
% point_cloud = point_cloud.pentagonsamples;
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
