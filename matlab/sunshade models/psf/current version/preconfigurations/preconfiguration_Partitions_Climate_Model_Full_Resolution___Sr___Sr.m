
function    [model_parameters] = ...
			........................................................................................................................................................
            preconfiguration_Partitions_Climate_Model_Full_Resolution___Sr___Sr (model_parameters)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Define parameters to analyze irradiance factors across the full sun-facing disc of the Earth, producing a 192-latitude x 288-longitude CESM-compatible   %&%%
%%%%	output. The disc spans the sunlit hemisphere (longitude ±90° from the sub-stellar point). The 144 columns on the nightside are set to 1.0 in              %&%%
%%%%    analysis_General___0___X. The sub-stellar point is placed at CESM longitude 0°; the CESM science team handles rotating the Earth underneath this         %&%%
%%%%    fixed sun-frame pattern when applying the file in a climate run.                                                                                         %&%%
%%%%                                                                                                                                                             %&%%
%%%%    Compared to 'preconfigured: low resolution', which computes only the central longitude strip (theta.n = 1) and replicates it across all 288              %&%%
%%%%    CESM longitudes, this configuration computes all 144 disc columns (theta.n = 144) so that each longitude bin receives its physically correct             %&%%
%%%%    irradiance factor based on its angular distance from the sub-stellar point.                                                                              %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the desired analysis type: 'irradiance', 'shaded irradiance', 'shading factor', or 'irradiance factor'.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


model_parameters.analysis_type               = 'irradiance factor';


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define partition half-steps corresponding to the output data format.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


half_step_latitude                           = 0.942408376963357/2;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the parameters of the surface partitions of the Star: angle intervals (in degrees), numbers of subintervals to partition.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


model_parameters.star.phi.interval           = [-90,90];
model_parameters.star.theta.interval         = [-90,90];

model_parameters.star.phi.n                  = 21;
model_parameters.star.theta.n                = 21;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the parameters of the surface partitions of the Planet.
%
%   Latitude (phi): 192 rows spanning -90° to +90°, with half-step extensions so that bin centres align with the CESM f09 latitude grid.
%
%   Longitude (theta): 144 columns spanning -90° to +90°.
%       - Bin centres fall at ±0.625°, ±1.875°, …, ±89.375° (1.25° spacing, matching the CESM longitude bin width).
%       - Bins 1–72   cover θ = -89.375° to -0.625°  (western dayside, mapped to CESM lon 270°–358.75° in analysis_General).
%       - Bins 73–144 cover θ = +0.625°  to +89.375° (eastern dayside, mapped to CESM lon   0°– 88.75° in analysis_General).
%
%   Disc validity: a cell at (phi, theta) lies inside the Earth's disc when sin²(phi) + sin²(theta) ≤ 1.
%   Corner cells that exceed this boundary (large |phi| and large |theta| simultaneously) are outside the
%   Earth's edge and the irradiance computation returns 1.0 for them naturally — no special handling needed.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


model_parameters.planet.phi.interval         = [-90 - half_step_latitude , 90 + half_step_latitude];
model_parameters.planet.theta.interval       = [-90, 90];

model_parameters.planet.phi.n                = 192;
model_parameters.planet.theta.n              = 144;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input strings indicating the type of Star and Planet surface partitions {'sphere', 'disc'}, and that of the Shade {'heliogyros'}.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


model_parameters.star.partition_type         = 'sphere';
model_parameters.planet.partition_type       = 'disc';
model_parameters.shade.partition_type        = 'heliogyros';


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define regular and leap year start times, frequencies, and number of periods.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


model_parameters.time.regular_year.start     = '2025-01-01T12:00:00';
model_parameters.time.regular_year.frequency = 1;
model_parameters.time.regular_year.periods   = 365;


model_parameters.time.leap_year.start        = '2024-01-01T12:00:00';
model_parameters.time.leap_year.frequency    = 1;
model_parameters.time.leap_year.periods      = 366;


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
