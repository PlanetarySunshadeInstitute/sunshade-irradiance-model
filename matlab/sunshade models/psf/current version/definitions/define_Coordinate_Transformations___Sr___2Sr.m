
function    [change_of_bases, origins] = ...
			........................................................................................................................................................ 
            define_Coordinate_Transformations___Sr___2Sr(kinematics_SP_data) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Define model elements with respect to a single arbitrary cartesian coordinate system                                                                    %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Represent the origins of the coordinate systems used for the Planet, Star, and Shade with respect to a single arbitrary coordinate system. Here we're choosing: 
%
%		- a standard Cartesian system
%		- origin at the Planet's center
%	 	- x axis that passes through the Planet and Star centers, increasing in the direction of the Star
%		- y axis in the orbital plane, and z axis the resulting normal
%
%	The unit vectors for the Star are then rotated 180 degrees about the z axis (with basis_S, which also represents the relevant change of basis matrix), and its
%	displacement is identified with its origin.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


change_of_bases.planet = [1 0 0 ; 0 1 0 ; 0 0 1];
change_of_bases.star   = [-1 0 0 ; 0 -1 0 ; 0 0 1];
change_of_bases.shade  = [1 0 0 ; 0 1 0 ; 0 0 1];

origins.planet         = [0 ; 0 ; 0];
origins.star           = [kinematics_SP_data.distance.SP ; 0 ; 0];
origins.shade          = [kinematics_SP_data.distance.PL1 ; 0 ; 0];


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%