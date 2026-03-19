
function    [partition_normal_vectors] = ...
			........................................................................................................................................................
            partition_spherical_unit_normal_vectors___N___M



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Define the lower and upper bounds for the standard spherical coordiantes (theta, phi) on the unit sphere, and the resolution of their partitions.       %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define lower and upper bounds for theta and phi, and calculate the magnitude of their resulting ranges.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


theta_lower_bound        = -90;
theta_upper_bound        = 90;
theta_range              = theta_upper_bound - theta_lower_bound;

phi_lower_bound          = -90;
phi_upper_bound          = 90;
phi_range                = phi_upper_bound - phi_lower_bound;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define the number of elements desired in each partition.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


n_theta                  = 8;
n_phi                    = 8; 


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Create an (n_theta x n_phi x 3) matrix containing the cartesian coordinates of the normal vectors at the center of each partitioned surface.            %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the partitions of theta and phi. To position each normal vector in the middle of its partitioned region, each angle partition stops a half-increment
%   after the lower bound, and finishes a half-increment before the upper bound (the latter being halted natively by the Matlab colon syntax).    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


increment_theta          = theta_range / n_theta;
partition_theta          = (theta_lower_bound + increment_theta / 2) : increment_theta : theta_upper_bound;

increment_phi            = phi_range / (n_phi);
partition_phi            = (phi_lower_bound + increment_phi / 2) : increment_phi : phi_upper_bound;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	For each partitioned (theta,phi), convert the normal vector to sphereical coordinates.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_normal_vectors = zeros(n_theta, n_phi, 3);


for i = 1:n_theta

	for j = 1:n_phi


		partition_normal_vectors(i,j,1:3) = [
						  					  sind(partition_phi(1,j))*sind(partition_theta(1,i)) ; 
						  					  sind(partition_phi(1,j))*cosd(partition_theta(1,i)) ;
						  					  cosd(partition_phi(1,j)) 
						 					];

	end

end


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%