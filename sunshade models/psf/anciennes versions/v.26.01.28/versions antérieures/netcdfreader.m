% plot_shading_t0.m
% Read instantaneous shading (t=0) from NetCDF and plot a world map.

clear; clc;

% --- File (assumed to be in current folder)
ncfile = 'instantaneous_shading_t0.nc';
assert(isfile(ncfile), 'File not found: %s', ncfile);

% --- Read data
lat   = double(ncread(ncfile, 'lat'));            % [-90..90], column vector
lon   = double(ncread(ncfile, 'lon'));            % [-180..179], row vector
shade = double(ncread(ncfile, 'shading_percent'));% size = [numel(lat) x numel(lon)]

% Optional: fetch a few global attributes for the title
tit = tryreadatt(ncfile, '/', 'title');
sumry = tryreadatt(ncfile, '/', 'summary');

% --- Stats
globalMean = mean(shade(~isnan(shade)), 'all');
cmax = max(shade(:), [], 'omitnan');  if isempty(cmax) || ~isfinite(cmax), cmax = 1; end
cmax = max(cmax, 1e-6); % avoid zero color range

% --- Plot
hasMapping = exist('axesm','file')==2 && license('test','map_toolbox');
figure('Color','w');

if hasMapping
    % Mapping Toolbox branch (Mercator, clipped to avoid singularities)
    latMaxPlot = 89;
    [LON, LAT] = meshgrid(lon, lat);
    LATclip = max(min(LAT, latMaxPlot - 1e-6), -latMaxPlot + 1e-6);

    axesm('mercator', ...
          'MapLatLimit',[-latMaxPlot latMaxPlot], ...
          'MapLonLimit',[-180 180], ...
          'Frame','on','Grid','on', ...
          'MeridianLabel','on','ParallelLabel','on', ...
          'PLineLocation',30,'MLineLocation',60, ...
          'FLineWidth',0.5);
    axis off; hold on;

    pcolorm(LATclip, LON, shade);
    shading flat; colormap(turbo); caxis([0 cmax]);
    cb = colorbar; cb.Label.String = '% sunlight reduction';

    % Coastlines if available
    try
        land = shaperead('landareas','UseGeoCoords',true);
        geoshow([land.Lat], [land.Lon], 'Color','k', 'LineWidth',0.8);
    catch
        % fall back silently if shapefile is not available
    end
else
    % Fallback: regular lat-lon image (Plate Carrée)
    imagesc(lon, lat, shade);
    set(gca,'YDir','normal');
    axis tight; box on;
    xlabel('Longitude (degrees)'); ylabel('Latitude (degrees)');
    colormap(turbo); caxis([0 cmax]);
    cb = colorbar; cb.Label.String = '% sunlight reduction';
    grid on; hold on;

    % Coastlines (if the 'coastlines' MAT file exists)
    if exist('coastlines.mat','file') || exist('coastlines','var')
        try
            S = load('coastlines'); % provides S.coastlat, S.coastlon
            plot(S.coastlon, S.coastlat, 'k', 'LineWidth', 0.5);
        catch
            % ignore if not available
        end
    end
end

% --- Title
titleStr = 'Instantaneous Sunlight Reduction (%), t=0';
if ~isempty(tit), titleStr = tit; end
subtitleStr = sprintf('Global mean = %.3g %%', globalMean);
if ~isempty(sumry), subtitleStr = sprintf('%s\n%s', sumry, subtitleStr); end
title({titleStr; subtitleStr});

% --- Save a PNG next to the NetCDF
outPng = 'instantaneous_shading_t0.png';
exportgraphics(gcf, outPng, 'Resolution', 200);
fprintf('Saved figure: %s\n', outPng);

% ---------- helper: safe attribute read ----------
function val = tryreadatt(ncfile, loc, name)
    val = '';
    try
        v = ncreadatt(ncfile, loc, name);
        if isstring(v) || ischar(v), val = char(v); end
    catch
        % attribute missing
    end
end
