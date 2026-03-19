
function    [position_and_velocity_Earth] = ...
			........................................................................................................................................................
            position_and_velocity_Earth___St___V (time_UTC)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Calculate position and velocity of Earth, at the input time, in a specified reference frame                                                             %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define relevant path to MuSCAT and load relevant SPICE kernels.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


MuSCAT_SPICE_add_paths_and_load_kernels___0___0


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Convert input UTC time (e.g. '2025-01-01T12:00:00') to "ephemeris time" (ET).
%
%	"CSPICE_STR2ET converts a string representing an epoch to a double precision value representing the number of TDB seconds past the J2000 epoch corresponding 
%	to the input epoch."
%
%   (https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/MATLAB/mice/cspice_str2et.html)
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


time_ET                     = cspice_str2et(time_UTC);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	"CSPICE_SPKEZR returns the state [position (km) and velocity (km/s)] of a target body relative to an observing body, optionally corrected for light time  
%	(planetary aberration) and stellar aberration." 
%
%	(https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/IDL/icy/cspice_spkezr.html)   
%
%	Inputs: 	(name of target body, ephemeris time, name of reference frame, aberration correction, name of observing body)
%	Outputs:	(6x1 vector indicating position [1:3,1] and velocity [4:6,1] in Cartesian coordinates)	
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


position_and_velocity_Earth = cspice_spkezr('EARTH', time_ET, 'J2000', 'NONE', 'SUN');


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%