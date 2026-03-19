
function    [solar_irradiances_P] = ...
			........................................................................................................................................................ 
            analysis_Solar_Irradiance_Spherical_Coordinate_Partitions___0___M 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Input: model partition and time parameters                                                                                                              %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Start program stop-watch.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


tic


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the parameters of the spherical coordinate partitions of the Star: angle intervals (in degrees), numbers of subintervals to partition.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


phi_interval_S   = [-90,90];
theta_interval_S = [-90,90];

n_phi_S          = 20;
n_theta_S        = 20;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the parameters of the spherical coordinate partitions of the Planet: angle intervals (in degrees), numbers of subintervals to partition.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


phi_interval_P   = [-90,90];
theta_interval_P = [-90,90];

n_phi_P          = 11;
n_theta_P        = 11;


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


[ ~ , radius_S, ~ , radius_P, ~ , luminosity_S] = constants_Physical_And_Estimated___0___R;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the relative Planet-Star position vector (in km) at the desired day from SPICE data saved in the Excel file, then calculate it's norm.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


relative_position_PS                            = import_Sun_Earth_Position___St___C (time_UTC);
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
%	Represent the basis and origin locations of the Planet and Star with respect to a single arbitrary coordinate system. Here we're choosing: 
%
%		- a standard Cartesian system
%		- origin at the Planet's center
%	 	- x axis that passes through the Planet and Star centers, increasing in the direction of the Star
%		- y axis in the orbital plane, and z axis normal to it
%
%	The unit vectors for the Star are then rotated 180 degrees about the z axis (with basis_S, which also represents the relevant change of basis matrix), and its
%	displacement is identified with its origin.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


basis_P                                           = [1 0 0 ; 0 1 0 ; 0 0 1];
basis_S                                           = [-1 0 0 ; 0 -1 0 ; 0 0 1];

origin_P                                          = [0 ; 0 ; 0];
origin_S                                          = [relative_position_norm_PS ; 0 ; 0];


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Calculate: intermediate and final results																			                                    %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Identify the differences between all partitioned surface points on the Star and Planet, and the projections of these onto the normal vectors in these sets. 
%	The relevant results are all scalars output in arrays whose first dimensions correspond to the Star and latter to the Planet: {phi_S, theta_S, phi_P, theta_P}. 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[norms_SP, projections_SPxS, projections_PSxP] = norms_BA_Projections_ABxA_BAxB___4M2C2S___2M...
									             (normal_vectors_S, normal_vectors_P, basis_S, basis_P, origin_S, origin_P, radius_S, radius_P);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Use one factor of the Star's projections to account for surface alignment, and another expanded in a power series to account for limb darkening.       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


projections_SPxS                               = projections_SPxS .* power_Series_Of_Matrix___MC___M (projections_SPxS, [0.3 ; 0.93 ; -0.23]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the surface area of each partitioned sub-surface on the Star (for use with output flux). By symmetry, this is only dependent on phi_S.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partitioned_area_ratios_S                      = partition_Area_Ratios_Unit_Sphere___2R2S___M(phi_interval_S, theta_interval_S, n_phi_S, n_theta_S);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the solar irradiance incident on each partitioned Planet sub-surface from each partitioned Star sub-surface.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


solar_irradiances_P                            = projections_SPxS * luminosity_S .* partitioned_area_ratios_S .* norms_SP.^(-2) .* projections_PSxP;
 

%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Identify the number of non-Matlab-trivialized input dimensions (>1) in the defined Planet partition so it's preserved when summing irradiance contributions.  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


non_trivial_dimensions_P                       = identify_Non_Trivial_Dimensions___2S___S (n_phi_P, n_theta_P);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Sum solar irradiance contributions over the Star's partition (all non-trival {phi_S, theta_S} dimensions) to get total solar irradiance at each element of the
%	planet partition (indexed by {phi_P, theta_P}). Then remove any trivial dimensions left by default of the "sum" function, and convert irradiance units from 
%   km^-2 to m^-2.       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


solar_irradiances_P                            = sum(solar_irradiances_P, [1 : ndims(solar_irradiances_P)-non_trivial_dimensions_P]);
solar_irradiances_P                            = squeeze(solar_irradiances_P);
solar_irradiances_P                            = solar_irradiances_P * (10^(-3))^2;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	End program stop-watch.  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


toc


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


