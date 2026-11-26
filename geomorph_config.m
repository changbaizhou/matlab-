function cfg = geomorph_config()
%GEOMORPH_CONFIG Reproducible geographic limits and metric analysis settings.
cfg.epsg=32648;
cfg.cellSize_m=30;
cfg.outlet_lonlat=[102.2287 37.9285]; % Reservoir dam vicinity; snapped to DEM network.
cfg.station_lonlat=[102.05 37.8666667]; % Independent Jiutiaoling area cross-check.
cfg.threshold_km2=1;
cfg.theta=0.45;
cfg.referenceArea_m2=1e6;
cfg.reachLength_m=1000;
cfg.reliefWindow_cells=101; % 3030 m square; computed before study-area masking.
cfg.excludeOutletReach_m=2500; % Avoid reservoir and dam influence in channel indices.
cfg.knickpointTolerance_m=30;
cfg.mountainBoundaryFile=fullfile('results','mountain_boundary_utm.csv');
cfg.swath_longitudes=[101.79 101.88 101.98];
cfg.swathWidth_m=3000;
cfg.swathBin_m=300;
cfg.thetaSensitivity=[0.35 0.45 0.55];
cfg.thresholdSensitivity_km2=[0.5 1 2];
end
