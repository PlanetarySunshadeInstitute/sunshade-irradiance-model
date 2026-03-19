function    [partition_areas_phi_dependent] = ...
			........................................................................................................................................................
            partition_Unit_Sphere_Into_Areas___VSS___M (phi_bounds_lower_and_upper, n_theta, n_phi)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Return a vector containing the areas of each sub-surface as a function of phi (by symmetry, all sub-surfaces with fixed phi will have the same area).   %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Partition the input intervals for theta and phi into n_theta and n_phi equal parts, respectively, then return vectors identifying the boundaries of each 
%	subinterval. 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_phi_boundaries               = partition_1D_Interval_To_N_Subinterval_Boundaries___VS___V(phi_bounds_lower_and_upper, n_phi);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define a vector to store partition variables as a function of the sphereical coordinate phi (as by symmetry, sub-surfaces with fixed phi would have the 
%	same area).      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_areas_phi_dependent          = zeros(1,n_phi);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the area of each sub-surface using:
%
%		- the result of the double integral in sphereical coordinates to get the total area of sub-surfaces with fixed phi. 
%		- symmetry to divide these into n_theta sub-surfaces of equal magnitude.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


for i = 1:n_phi

	partition_areas_phi_dependent(1,i) = 2*pi*abs(sind(partition_phi_boundaries(1,i+1)) - sind(partition_phi_boundaries(1,i)));

end


partition_areas_phi_dependent          = partition_areas_phi_dependent ./ n_theta;


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%