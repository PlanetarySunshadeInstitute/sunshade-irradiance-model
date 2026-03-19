
function    [unit_normal_vectors] = ...
			........................................................................................................................................................
            partition_Unit_Normal_Vectors_On_Unit_Sphere___2R2S___M (phi_interval, theta_interval, n_phi, n_theta)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Create an (3 x n_phi x n_theta) matrix containing the cartesian coordinates of the normal vectors at the center of each partitioned sub-surface.        %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Partition the input intervals for phi and theta into n_phi and n_theta equal parts, respectively, then identify the centers of each subinterval. 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_phi_centers                = generate_1D_Subinterval_Centers___RS___R(phi_interval, n_phi);
partition_theta_centers              = generate_1D_Subinterval_Centers___RS___R(theta_interval, n_theta);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the unit normal vector (in cartesian coordinatse) at the center of each subsurface. Save vector information along the resulting arrays 1st 
%	dimension, and use its 2nd and 3rd dimensions to index phi and theta, respectively.  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


unit_normal_vectors                  = zeros(3, n_phi, n_theta);


for i = 1:n_phi

	for j = 1:n_theta


		unit_normal_vectors(1:3,i,j) = [
						  			    cosd(partition_phi_centers(1,i))*cosd(partition_theta_centers(1,j)) ; 
						  			    cosd(partition_phi_centers(1,i))*sind(partition_theta_centers(1,j)) ;
						  			    sind(partition_phi_centers(1,i)) 
						 			   ];

	end

end


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%