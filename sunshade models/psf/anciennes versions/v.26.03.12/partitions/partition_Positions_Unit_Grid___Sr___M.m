
function    [scaling_factors] = ...
			........................................................................................................................................................
            partition_Positions_Unit_Grid___Sr___M (partition_parameters)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Create a 3xn matrix containing disc scaling factors (that when multiplied on the right by the disc's normal, will produce position vectors).            %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Partition the input intervals for phi and theta into n_phi and n_theta equal parts, respectively, then identify the centers of each subinterval. 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_phi_centers            = generate_1D_Subinterval_Centers___RS___R(partition_parameters.phi.interval, partition_parameters.phi.n);
partition_theta_centers          = generate_1D_Subinterval_Centers___RS___R(partition_parameters.theta.interval, partition_parameters.theta.n);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the midpoint positions of 3D latitude/longitude distances from their respective center axis.  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


scaling_factors                  = zeros(3, partition_parameters.phi.n, partition_parameters.theta.n);


for i = 1:partition_parameters.phi.n

	for j = 1:partition_parameters.theta.n


		scaling_factors(1:3,i,j) = [
						  			    0 ; 
						  			    sind(partition_theta_centers(1,j)) ;
						  			    sind(partition_phi_centers(1,i)) ;
						 			   ];

	end

end


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%