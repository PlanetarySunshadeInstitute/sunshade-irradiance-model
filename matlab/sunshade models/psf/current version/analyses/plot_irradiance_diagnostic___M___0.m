
function plot_irradiance_diagnostic___M___0 (results_2d, lat_deg, time_doy, save_path, input_filename, num_shades)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Plot a two-panel irradiance factor diagnostic figure:                                                                                                   %&%%
%%%%	  Panel 1 (left)  — degrees from sub-stellar point vs month contour chart.                                                                            %&%%
%%%%	  Panel 2 (right) — scalar summary: constellation metadata and area-weighted                                                                           %&%%
%%%%	                    mean irradiance factor for each of the 12 diagnostic months.                                                                       %&%%
%%%%                                                                                                                                                            %&%%
%%%%	Both panels share the same figure window. In preconfigured mode the entire                                                                             %&%%
%%%%	figure is saved as a single PNG alongside the NC output file.                                                                                          %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   INPUTS
%     results_2d      [N_lat x N_time]   irradiance factor values (central longitude column, monthly samples)
%     lat_deg         [1 x N_lat]        degrees from sub-stellar point, e.g. linspace(-90, 90, N_lat)
%     time_doy        [1 x N_time]       day-of-year for each time column (1-indexed)
%     save_path       (optional string)  full path incl. filename for PNG output;
%                                        if omitted the figure is left open for inspection
%     input_filename  (optional string)  name of the constellation .mat input file
%     num_shades      (optional scalar)  number of shade objects in the constellation
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if nargin < 4, save_path      = '';  end
if nargin < 5, input_filename = '';  end
if nargin < 6, num_shades     = [];  end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Month tick positions and labels.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


month_tick_doy = [1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335];
month_labels   = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Colormap: blue (fully shaded) → white (no shading).  Matches S&M convention.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


n_colors = 256;
frac     = linspace(0, 1, n_colors)';
r_ch     = interp1([0, 1], [0.10, 1.00], frac);
g_ch     = interp1([0, 1], [0.30, 1.00], frac);
b_ch     = interp1([0, 1], [0.80, 1.00], frac);
cmap     = [r_ch, g_ch, b_ch];


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Contour level range.  Guard against flat (all-identical) input.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


v_min = min(results_2d(:));
v_max = max(results_2d(:));
if abs(v_max - v_min) < 1e-6
	v_max = v_min + 0.01;
end

n_levels = 16;
levels   = linspace(v_min, v_max, n_levels + 1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Wrap one column at each edge so the chart fills Jan-to-Dec without blank margins.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


time_ext = [time_doy(end) - 365,  time_doy(:)',  time_doy(1) + 365];
data_ext = [results_2d(:, end),   results_2d,    results_2d(:, 1)  ];

[T, L] = meshgrid(time_ext, lat_deg(:)');


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Scalar metrics — area-weighted monthly means and overall mean.
%
%	  weights        cos-latitude area weights (sum to 1).
%	  monthly_means  area-weighted mean irradiance factor for each time column.
%	  overall_mean   simple mean of the monthly area-weighted means.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


weights       = cosd(lat_deg(:));
weights       = weights / sum(weights);
monthly_means = weights' * results_2d;       % 1 × N_time
overall_mean  = mean(monthly_means);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Build the summary text for Panel 2.
%
%	Constellation block: filename and shade count (shown as-is if available,
%	otherwise a dash placeholder).
%
%	Monthly block: one row per time column, labelled by the nearest month name.
%	A separator line precedes the overall mean.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if ~isempty(input_filename)
	filename_str = input_filename;
else
	filename_str = '—';
end

if ~isempty(num_shades)
	shades_str = sprintf('%d shades', num_shades);
else
	shades_str = '—';
end

monthly_rows = '';
for k = 1 : numel(time_doy)
	[~, m_idx]   = min(abs(month_tick_doy - time_doy(k)));
	monthly_rows = [monthly_rows, ...
	                sprintf('  %-4s  %.6f\n', month_labels{m_idx}, monthly_means(k))]; %#ok<AGROW>
end

summary_str = sprintf([ ...
	'Constellation\n',                         ...
	'  %s\n',                                  ...
	'  %s\n',                                  ...
	'\n',                                      ...
	'Area-weighted mean DF\n',                 ...
	'%s',                                      ...
	'  ------------------\n',                  ...
	'  Mean  %.6f\n',                          ...
	], filename_str, shades_str, monthly_rows, overall_mean);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Create figure and two-panel tiled layout.
%	Both panels are approximately square and sit side by side.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


fig = figure ( ...
	'Name',        'Dimming Factor — Diagnostics', ...
	'NumberTitle', 'off', ...
	'Color',       'white', ...
	'Position',    [100, 100, 1060, 460] );

tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Panel 1: Latitude vs Month.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


ax1 = nexttile(tl);
hold(ax1, 'on');

[~, h_fill] = contourf(ax1, T, L, data_ext, levels);
set(h_fill, 'LineColor', [0.35 0.35 0.35], 'LineWidth', 0.4);

colormap(ax1, cmap);
clim(ax1, [v_min, v_max]);

cb1              = colorbar(ax1, 'eastoutside');
cb1.Label.String = 'Dimming factor';
cb1.Label.FontSize = 10;
cb1.FontSize       = 9;
cb1.TickLabelsMode = 'manual';
cb1.TickLabels     = arrayfun(@(v) sprintf('%.4f', v), cb1.Ticks, 'UniformOutput', false);

ylabel(ax1, 'Degrees from sub-stellar point [deg]', 'FontSize', 11);
title(ax1,  'Dimming factor — degrees from sub-stellar point vs month', 'FontSize', 12, 'FontWeight', 'bold');

xlim(ax1, [1, 365]);
ylim(ax1, [-90, 90]);
yticks(ax1,      [-90, -60, -30, 0, 30, 60, 90]);
yticklabels(ax1, {'-90°', '-60°', '-30°', '0°', '30°', '60°', '90°'});
xticks(ax1,      month_tick_doy);
xticklabels(ax1, month_labels);
set(ax1, 'TickDir', 'out', 'Box', 'on', 'FontSize', 10, 'Layer', 'top', ...
         'XGrid', 'off', 'YGrid', 'off');

hold(ax1, 'off');


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Panel 2: Summary statistics.
%
%	The panel is a plain text box with no visible axes.  Metrics are derived from
%	the same 12-monthly central-column sample used by Panel 1.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


ax2 = nexttile(tl);
axis(ax2, 'off');
title(ax2, 'Summary statistics', 'FontSize', 12, 'FontWeight', 'bold', 'Visible', 'on');

text(ax2, 0.05, 0.88, summary_str, ...
	'Units',               'normalized', ...
	'VerticalAlignment',   'top', ...
	'HorizontalAlignment', 'left', ...
	'FontSize',            10, ...
	'FontName',            'Courier New', ...
	'Interpreter',         'none');


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Render.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


drawnow;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Save to PNG (preconfigured mode) or leave open for inspection (manual mode).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if ~isempty(save_path)
	figure(fig);
	drawnow;
	try
		exportgraphics(fig, save_path, 'Resolution', 300);
	catch
		try
			saveas(fig, save_path, 'png');
		catch err
			warning('%s', sprintf('Could not save diagnostic chart: %s\n  %s', save_path, err.message));
		end
	end
	fprintf('Saved diagnostic chart: %s\n', save_path);
	close(fig);
end


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
