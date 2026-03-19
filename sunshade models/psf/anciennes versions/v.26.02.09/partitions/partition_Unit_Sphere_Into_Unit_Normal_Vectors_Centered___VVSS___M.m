
function    [partition_unit_normal_vectors_centered] = ...
			........................................................................................................................................................
            partition_Unit_Sphere_Into_Unit_Normal_Vectors_Centered___VVSS___M (theta_bounds_lower_and_upper, phi_bounds_lower_and_upper, n_theta, n_phi)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Create an (n_theta x n_phi x 3) matrix containing the cartesian coordinates of the normal vectors at the center of each partitioned sub-surface.        %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Partition the input intervals for theta and phi into n_theta and n_phi equal parts, respectively, then return vectors identifying the centers of each 
%	subinterval. 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_theta_centers                                 = partition_1D_Interval_To_N_Subinterval_Centers___VS___V(theta_bounds_lower_and_upper, n_theta);
partition_phi_centers                                   = partition_1D_Interval_To_N_Subinterval_Centers___VS___V(phi_bounds_lower_and_upper, n_phi);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	For (theta,phi) at the center of each partitioned sub-surface, calculate the unit normal vector in Cartesian coordinates.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_unit_normal_vectors_centered                  = zeros(n_theta, n_phi, 3);


for i = 1:n_theta

	for j = 1:n_phi


		partition_unit_normal_vectors_centered(i,j,1:3) = [
						  					               cosd(partition_phi_centers(1,j))*sind(partition_theta_centers(1,i)) ; 
						  					               cosd(partition_phi_centers(1,j))*cosd(partition_theta_centers(1,i)) ;
						  					               sind(partition_phi_centers(1,j)) 
						 					              ];

	end

end


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%