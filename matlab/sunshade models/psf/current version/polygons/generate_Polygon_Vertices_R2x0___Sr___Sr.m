
function    [polygon] = ...
			........................................................................................................................................................
            generate_Polygon_Vertices_R2x0___Sr___Sr (parameters) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Generate the input regular polygon as a 2D polyshape in Matlab.                                                                                         %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Partition the inner and outer angles (2pi) based on the input number of vertices.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


angle_partitions_inner           = [0 : (2*pi/parameters.sides.inner) : 2*pi];
angle_partitions_outer           = [0 : (2*pi/parameters.sides.outer) : 2*pi];


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	If the partition created a point at 2pi, remove this point (as it's covered by 0).     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if angle_partitions_inner(1,end) == 2*pi

	angle_partitions_inner       = angle_partitions_inner(1,1:end-1);

end


if angle_partitions_outer(1,end) == 2*pi

	angle_partitions_outer       = angle_partitions_outer(1,1:end-1);

end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Generate the sets of coordinates describing the inner and outer contours.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


polygon.R2x0.vertices.inner(1,:) = parameters.radii.inner * cos(angle_partitions_inner);
polygon.R2x0.vertices.inner(2,:) = parameters.radii.inner * sin(angle_partitions_inner);
polygon.R2x0.vertices.inner(3,:) = 0;


polygon.R2x0.vertices.outer(1,:) = parameters.radii.outer * cos(angle_partitions_outer);
polygon.R2x0.vertices.outer(2,:) = parameters.radii.outer * sin(angle_partitions_outer);
polygon.R2x0.vertices.outer(3,:) = 0;


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%