function    [partition_area_ratios] = ...
			........................................................................................................................................................
            partition_Area_Ratios_Unit_Sphere___2R2S___M (phi_interval, theta_interval, n_phi, n_theta)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Return a vector containing partition area ratios (sub-surface : unit sphere) as a function of phi (by symmetry, these are independent of theta).        %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Partition the input intervals for theta and phi into n_theta and n_phi equal parts, respectively, and return a vector identifying the boundaries of each 
%	subinterval. Then create vectors identifying the upper and lower boundaries.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_phi_boundaries       = generate_1D_Subinterval_Boundaries___VS___V(phi_interval, n_phi);
partition_phi_boundaries_upper = partition_phi_boundaries(2:end);
partition_phi_boundaries_lower = partition_phi_boundaries(1:end-1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Take the ratio of the interval for theta over its full 360 degree range.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


theta_interval_ratio           = abs(theta_interval(1,1) - theta_interval(1,2)) / 360;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the partition sub-surface area ratios to that of the unit sphere by: 
%
%		- finding the areas of the strips [phi_i.phi_(i+1)]x[0,360] via the double integral in sphereical coordinates (with phi giving the angle to the xy plane)
%       - changing the coefficient in this integral from 2pi to 2pi/4pi = 1/2 to represent the ratio of the strips' areas to that of the unit sphere
%		- (by symmetry) multiply each strip by the ratio (theta_interval):360 to get the area of the union of partion sub-surfaces within each stip  
%		- (by symmetry) divide these ratios into n_theta sub-surfaces of equal magnitude.  
%
%		- convert the results to a column vector so that phi varies along the first dimension.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_area_ratios          = (1/2) * abs(sind(partition_phi_boundaries_upper) - sind(partition_phi_boundaries_lower)) * theta_interval_ratio / n_theta;
partition_area_ratios          = partition_area_ratios';


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%