
function    [solar_irradiances_P] = ...
			........................................................................................................................................................ 
            analysis_Solar_Irradiance_Spherical_Coordinate_Partitions___0___M 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Input: model partition and time parameters                                                                                                              %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the parameters of the spherical coordinate partitions of the Star: angle intervals (in degrees), numbers of subintervals to partition.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


phi_interval_S   = [-90,90];
theta_interval_S = [-90,90];

n_phi_S          = 10;
n_theta_S        = 10;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the parameters of the spherical coordinate partitions of the Planet: angle intervals (in degrees), numbers of subintervals to partition.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


phi_interval_P   = [-90,90];
theta_interval_P = [-90,90];

n_phi_P          = 1;
n_theta_P        = 1;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the desired time in UTC     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


time_UTC         = '2025-01-01T12:00:00';


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Import: the radii of the Plant and Star (JPL, cited), the star's luminosity (cited), and the relative position of the Planet and Star (SPICE).          %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the estimated radii (km) of the Planet and the Star, and the star's luminosity (W). Calculate the star's luminosity per unit area (W/km^2).  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[ ~ , radius_S, ~ , radius_P, ~ , luminosity_S] = constants_Physical_And_Estimated___0___V;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import a relative Planet-Star position and velocity vector (in km, km/s), at the desired time, from SPICE. Then calculate the relative position's norm.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


relative_position_and_velocity_PS               = import_SPICE_Earth_Position_And_Velocity___St___V (time_UTC);
relative_position_PS                            = relative_position_and_velocity_PS(1:3,1);
relative_position_norm_PS                       = norm(relative_position_PS);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Generate: Model partitions, coordinate systems, and related vectors.                                                                                    %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Generate unit normal vectors centered on the sub-surfaces of the input unit sphere partitions to represent the Planet and the Star:
%		
%		- vectors are represented in standard cartesian coordiantes
%		- the first dimension stores vector information, the remaining dimensions index the partition (phi, theta)
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


normal_vectors_S                                  = partition_Unit_Normal_Vectors_On_Unit_Sphere___2R2S___M(phi_interval_S,theta_interval_S,n_phi_S,n_theta_S);
normal_vectors_P                                  = partition_Unit_Normal_Vectors_On_Unit_Sphere___2R2S___M(phi_interval_P,theta_interval_P,n_phi_P,n_theta_P);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Represent the basis and origin locations of the Planet and Star with respect to a single arbitrary coordinate system (and define its identity transformation). 
%	Here we're choosing:
%
%		- a standard Cartesian system
%		- origin at the planet's center
%	 	- x axis that passes through the Planet and Star centers, increasing in the direction of the star
%		- y axis in the orbital plane, and z axis normal to it
%
%	The unit vectors for the star are then rotated 180 degrees about the z axis (with basis_S, which also represents the relevant change of basis matrix), and its
%	displacement is identified with its origin.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


basis_P                                           = [1 0 0 ; 0 1 0 ; 0 0 1];
basis_S                                           = [-1 0 0 ; 0 -1 0 ; 0 0 1];
basis_I                                           = [1 0 0 ; 0 1 0 ; 0 0 1];

origin_P                                          = [0 ; 0 ; 0];
origin_S                                          = [relative_position_norm_PS ; 0 ; 0];


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Calculate: intermediate and final results																			                                    %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%
%%%%	Perform orientation-dependent calculations (independent of position).
%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Perform the linear transformation (projection via negative* dot product) of all pairs of Star and Planet normal vectors, and compose this with the indicator 
%	function on the positive reals*. *These two steps ensure flux is positive for surfaces facing each other, and zero otherwise. The scalar results are output 
%	whose first indices correspond to that of the Star's vectors, and last to that of the Planet's: (phi_S, theta_S, phi_P, theta_P).   
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_SxP                  = dot_products_In_AxB___4M___M (-normal_vectors_S, normal_vectors_P, basis_S, basis_P);
dot_products_in_SxP                  = max(0, dot_products_in_SxP);


%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%
%%%%	Perform calculations that depend on oritentation and position.
%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Take the differences between all vectors in normal_vectors_S and normal_vectors_P, multiplying unit normals by respective radii to get actual surface positions,
%	and accounting for differences in the origins of their representations. Then calculate the norms of these differences (last two arguments: take the L2 norm 
%	along the 1st dimension), and normalize the vectors.  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[norms_in_SxP, dot_products_in_SPxS] = norms_In_AxB_And_Dot_Products_In_ABxA___4M2C2S___2M...
									   (normal_vectors_S, normal_vectors_P, basis_S, basis_P, origin_S, origin_P, radius_S, radius_P);


%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%
%%%%	Calculate the area ratios of each partitioned sub-surface on the Star (ratios of total output flux). 
%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the surface area of each partitioned sub-surface on the Star (for use with output flux). By symmetry, this is only dependent on phi_S.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partitioned_area__ratios_S           = partition_Area_Ratios_Unit_Sphere___2R2S___M(phi_interval_S, theta_interval_S, n_phi_S, n_theta_S);


%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%
%%%%	Calculate the solar irradiance incident on each partitioned Planet sub-surface.
%%%%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the solar irradiance incident on each partitioned Planet sub-surface from each partitioned Star sub-surface.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


solar_irradiances_P                  = luminosity_S * dot_products_in_SxP .* partitioned_area__ratios_S .* norms_in_SxP.^(-2) .* dot_products_in_SPxS;
 

%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Identify the number of non-trivial dimensions associated with the Planet (so that they're preserved is summing over radiance contributions).  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


non_trivial_dimensions_P             = identify_Non_Trivial_Dimensions___2S___S (n_phi_P, n_theta_P);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Sum solar irradiance contributions over and non-trival {phi_S, theta_S} dimensions to et total solar irradiance indexed by {phi_P, theta_P} (each partitioned 
%   Planet sub-surface). Then remove any default trivial dimensions output from the "sum" function, and convert irradiance units from from km^-2 to m^-2.       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


solar_irradiances_P                  = sum(solar_irradiances_P, [1 : ndims(solar_irradiances_P)-non_trivial_dimensions_P]);
solar_irradiances_P                  = squeeze(solar_irradiances_P);
solar_irradiances_P                  = solar_irradiances_P * (10^(-3))^2;


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%