
function    [kinematics_SP_data] = ...
			........................................................................................................................................................
            distance_PL1___Sr___Sr(kinematics_SP_data)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	<Description>                                                                                                                                           %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the gravitational constant (N⋅km^2/kg^2) and the masses of the star and the planet (kg).     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[~ , masses , constant_gravitational, ~] = parameters_Physical___0___2Sr2S;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the angular velocity of the planet.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


angular_velocity_norm_P                  = angular_Velocity_Norm___2C___S (kinematics_SP_data.position, kinematics_SP_data.velocity);
%angular_velocity_norm_P                 = constant_gravitational*(masses.star+masses.planet)/kinematics_SP_data.distance.SP^3;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define the constant R and the star-planet distance and the L1 equilibrium equation. Then find the equations solution with an initial guess at the zero point.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


R                                        = kinematics_SP_data.distance.SP;
L1_equilibrium_equation                  = @(r) constant_gravitational*masses.star/(R-r)^2 - constant_gravitational*masses.planet/r^2 - angular_velocity_norm_P^2*(R-r);
kinematics_SP_data.distance.PL1          = fzero(L1_equilibrium_equation, 0.01 * kinematics_SP_data.distance.SP);


kinematics_SP_data.distance.PL1          = 2.36*10^6;
 

% L1_equilibrium = @(x) constant_gravitational*masses.star/(kinematics_SP_data.distance.SP-x)^2 - constant_gravitational*masses.planet/x^2 - (masses.star/(masses.star+masses.planet)*kinematics_SP_data.distance.SP - x)*(masses.star+masses.planet)/kinematics_SP_data.distance.SP^3;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	For reference, perform a simpler L1 calculation ignoring centripetal force.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


a                                        = masses.star - masses.planet;
b                                        = 2 * masses.planet * kinematics_SP_data.distance.SP;
c                                        = - masses.planet * kinematics_SP_data.distance.SP^2;


quadratic_estimate                       = (-b + sqrt(b^2-4*a*c))/(2*a);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%