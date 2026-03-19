
function    [kinematics_SP_data] = ...
			........................................................................................................................................................
            distance_PL1___Sr___Sr(kinematics_SP_data)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Import physical parameters and calculate the planet's angular velocity from the star planet kinematics data.                                            %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the gravitational constant (N⋅km^2/kg^2) and the masses of the star and the planet (kg).     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[~ , masses , constant_gravitational, ~] = parameters_Physical___0___2Sr2S;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the angular velocity of the planet using the Earth position (km) and velocity (km/s) vectors [result is in radians/s].      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


angular_velocity_norm_P                  = angular_Velocity_Norm___2C___S (kinematics_SP_data.position, kinematics_SP_data.velocity);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Solve for the distance bewtween the planet and the L1 point.                                                                                            %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define the constant R and the star-planet distance and the L1 equilibrium equation. Then find the equations solution with an initial guess at the zero point. 	
%   
%	Units: 
%
%			gravitational Constant: N⋅km^2/kg^2 = m⋅ km^2/kgs^2
%			masses: kg
%			distances: km
%			angular velocity: 1/s 
%
%	Note: the distance factor in the angular velocity term is converted from km to m to correspond with the uncancelled units (m) in the gravitational constant. 	  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


R                               = kinematics_SP_data.distance.SP;
L1_equilibrium_equation         = @(r) constant_gravitational*masses.star/(R-r)^2 - constant_gravitational*masses.planet/r^2 - angular_velocity_norm_P^2*(R-r)*10^3;


% ORIGINAL CODE: kinematics_SP_data.distance.PL1 = fzero(L1_equilibrium_equation, 0.01 * R); 
kinematics_SP_data.distance.PL1 = kinematics_SP_data.distance.SP * (masses.planet / (3 * masses.star))^(1/3);

%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%