
function    [model_parameters] = ...
			........................................................................................................................................................
            preconfiguration_Partitions_Climate_Model_Low_Resolution___Sr___Sr (model_parameters)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Define parameters to analyze irradiance factors according to a low resolution climate model output.                                                     %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the desired analysis type: 'irradiance', 'shaded irradiance', 'shading factor', or 'irradiance factor'.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


model_parameters.analysis_type               = 'irradiance factor';


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define partition half-steps corresponding to the output data format.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


half_step_latitude                           = 0.942408376963357/2;
half_step_longitude                          = 0.625;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the parameters of the surface partitions of the Star: angle intervals (in degrees), numbers of subintervals to partition.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


model_parameters.star.phi.interval           = [-90,90];
model_parameters.star.theta.interval         = [-90,90];

model_parameters.star.phi.n                  = 21;
model_parameters.star.theta.n                = 21;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%		Input the parameters of the surface partitions of the Planet: angle intervals (in degrees), numbers of subintervals to partition.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


model_parameters.planet.phi.interval         = [-90 - half_step_latitude , 90 + half_step_latitude];
model_parameters.planet.theta.interval       = [-half_step_longitude , half_step_longitude];

model_parameters.planet.phi.n                = 192;
model_parameters.planet.theta.n              = 1;


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
model_parameters.time.leap_year.periods      = 364;


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%