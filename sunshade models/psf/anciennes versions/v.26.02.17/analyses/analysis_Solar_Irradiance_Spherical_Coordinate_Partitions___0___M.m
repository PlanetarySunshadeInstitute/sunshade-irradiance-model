
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

n_phi_S          = 8;
n_theta_S        = 8;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Input the parameters of the spherical coordinate partitions of the Planet: angle intervals (in degrees), numbers of subintervals to partition.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


phi_interval_P   = [-90,90];
theta_interval_P = [-90,90];

n_phi_P          = 8;
n_theta_P        = 8;


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


normal_vectors_S                                  = partition_Unit_Normal_Vectors_On_Unit_Sphere___VVSS___M(phi_interval_S,theta_interval_S,n_phi_S,n_theta_S);
normal_vectors_P                                  = partition_Unit_Normal_Vectors_On_Unit_Sphere___VVSS___M(phi_interval_P,theta_interval_P,n_phi_P,n_theta_P);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Find the dimensions of the vectors contained within each set (these must be equal).      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


normal_vector_dimensions                          = size(normal_vectors_S, 1);


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


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	- Record the magnitudes of the original indice dimensions of the sets of normal vectors.
%	- Convert these sets to a 2D array (such that the mapping of indices in Z^n to Z cycles progressively through lower dimensional indices to higher: 
%	- e.g. (i,j) -> i+(j-1)|i|
%	- Perform a change of coordinates to the same coordinate system.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[normal_vectors_S, indice_dimension_magnitudes_S] = represent_Vectors_As_2D_Matrix (normal_vectors_S, basis_S);
[normal_vectors_P, indice_dimension_magnitudes_P] = represent_Vectors_As_2D_Matrix (normal_vectors_P, basis_P);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	<Description>     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


output_dimensions_SP                              = [indice_dimension_magnitudes_S indice_dimension_magnitudes_P];


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Perform orientation-dependent calculations (independent of position).                                                                           	    %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Perform the linear transformation (projection via negative* dot product) of all pairs of Star and Planet normal vectors, such that the first dimension of the 
%   resulting 2D matrix corresponds to the number of Star normal vectors (dimensions 2:end of normal_vectors_S) and the second to that of the planet (dimensions
%	2:end of normal_vectors_P). (*The negative sign ensures flux is positive when surface normal vectors face one another).   
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_SxP = -normal_vectors_S' * normal_vectors_P;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Reshape the array of (scalar results) so that the first dimensions correspond to the Star and the remaining to the planet: (phi_S, theta_S, phi_P, theta_P).     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_SxP = reshape(dot_products_in_SxP, output_dimensions_SP);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Compose this linear tranformation with the indicator function on the positive reals (to eliminate negative flux for surfaces not facing one another).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_SxP = max(0, dot_products_in_SxP);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Perform calculations that depend on oritentation and position.                                                                                          %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Permute the dimensions of the 2D array of the Planet's normal vectors so that they're now listed along its (once trivial) 3rd dimension. Then replicate the
%   resulting 2D (dimensions 1,3) sub-array along along the 2nd dimension to correspond to the magnitude of the 2nd dimensions of (number of vectors in) the star's
%	normal vectors.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


normal_vectors_P    = permute(normal_vectors_P,[1 3 2]);
normal_vectors_P    = repmat(normal_vectors_P, [1 length(normal_vectors_S(1,:)) 1]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Take the differences between all vectors in normal_vectors_S and normal_vectors_P, multiplying unit normals by respective radii to get actual surface positions,
%	and accounting for differences in the origins of their representations. Then calculate the norms of these differences (last two arguments: take the L2 norm 
%	along the 1st dimension), and normalize the vectors.  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


differences_in_SxP  = radius_S * normal_vectors_S - ((origin_P - origin_S) + radius_P * normal_vectors_P);

clear normal_vectors_P ;

norms_in_SxP        = vecnorm(differences_in_SxP, 2, 1);
differences_in_SxP  = differences_in_SxP ./ norms_in_SxP;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Normalize the difference vectors, then reshape their matrix representation to a 2D array. 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_SPxS = -differences_in_SxP .* normal_vectors_S;

clear differences_in_SxP normal_vectors_S


dot_products_in_SPxS = sum(dot_products_in_SPxS, 1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	reshape the array of (scalar results) so that the first dimensions correspond to the Star and the remaining to the planet: (phi_S, theta_S, phi_P, theta_P).      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_SPxS = reshape(dot_products_in_SPxS, output_dimensions_SP);
norms_in_SxP         = reshape(norms_in_SxP, output_dimensions_SP);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Clear all extraneous normal vector information.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Calculate the areas of each partitioned sub-surface on the Star, and the resulting flux hitting each partitioned sub-surface on the Planet.             %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the surface area of each partitioned sub-surface on the Star (for use with output flux). By symmetry, this is only dependent on phi_S.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partitioned_area__ratios_S   = partition_Area_Ratios_Unit_Sphere___VSS___M(phi_interval_S, theta_interval_S, n_phi_S, n_theta_S);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the solar irradiance incident on each partitioned Planet sub-surface from each partitioned Star sub-surface.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


solar_irradiances_P          = luminosity_S * dot_products_in_SxP .* partitioned_area__ratios_S .* norms_in_SxP.^(-2) .* dot_products_in_SPxS;
 

%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Identify the number of non-trivial dimensions associated with the Planet (so that they're preserved is summing over radiance contributions).  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


non_trivial_dimensions_P     = 0;


if n_phi_P > 1

	non_trivial_dimensions_P = non_trivial_dimensions_P + 1;

end


if n_theta_P > 1

	non_trivial_dimensions_P = non_trivial_dimensions_P + 1;

end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Sum solar irradiance contributions over and non-trival {phi_S, theta_S} dimensions to et total solar irradiance indexed by {phi_P, theta_P} (each partitioned 
%   Planet sub-surface). Then remove any default trivial dimensions output from the "sum" function, and convert irradiance units from from km^-2 to m^-2.       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


solar_irradiances_P          = sum(solar_irradiances_P, [1 : ndims(solar_irradiances_P)-non_trivial_dimensions_P]);
solar_irradiances_P          = squeeze(solar_irradiances_P);
solar_irradiances_P          = solar_irradiances_P * (10^(-3))^2;


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%