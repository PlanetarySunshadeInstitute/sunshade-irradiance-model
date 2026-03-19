
function    [partitions] = ...
			........................................................................................................................................................ 
            generate_Partitions___A___Sr(model_parameters, partition_bodies) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Import the input model parameters and generate the model partitions accordingly.                                                                        %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	If the Star is listed amongst the requested partitioned bodie, identify the Star's partition type and generate its unit normal vectors accordingly.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if any(ismember(partition_bodies,'star'))


	switch model_parameters.star.partition_type 

		case 'sphere'

			partitions.unit_normal_vectors.star   = partition_Unit_Normal_Vectors_On_Unit_Sphere___Sr___M(model_parameters.star);

	end


end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	If the Planet is listed amongst the requested partitioned bodie, identify the Star's partition type and generate its unit normal vectors accordingly.   
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if any(ismember(partition_bodies,'planet'))


	switch model_parameters.planet.partition_type 

		case 'sphere'

			partitions.unit_normal_vectors.planet = partition_Unit_Normal_Vectors_On_Unit_Sphere___Sr___M(model_parameters.planet);

	end


end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	If the Shade is listed amongst the requested partitioned bodie, identify the Star's partition type and generate its unit normal vectors accordingly.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if any(ismember(partition_bodies,'shade'))


	switch model_parameters.shade.partition_type 

		case 'heliogyros'

			partitions.shade.heliogyros           = import_Heliogyro_Kinematics_Data___0___Sr;

	end


end


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%