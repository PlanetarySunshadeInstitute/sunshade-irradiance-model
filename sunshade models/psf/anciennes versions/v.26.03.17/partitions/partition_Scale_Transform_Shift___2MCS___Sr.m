
function    [partition_position_vectors] = ...
			........................................................................................................................................................
            partition_Scale_Transform_Shift___2MCS___Sr (partition_position_vectors, change_of_basis, origin)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Perform a change of basis on the position and normal vectors, relocate the origin of the position vectors, and normalize the normal vectors.            %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Generate the identity matrix corresponding to the dimensions of the input vectors (the magnitude of the first dimension of their input matrix), record this
%	matrix's initial dimensions, then reshape it into a 2D matrix (where vector information is stored along the first dimension).      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


identity_transformation                                   = eye(size(partition_position_vectors,1));
[partition_position_vectors, indice_dimension_magnitudes] = represent_Vectors_As_2D_Matrix___2M___MR (partition_position_vectors, identity_transformation);
vector_dimensions                                         = size(partition_position_vectors,1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Perform the input change of basis transformation.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_position_vectors                                = change_of_basis * partition_position_vectors;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Shift the origin of the partition's coordinate system to that of the common.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_position_vectors                                = partition_position_vectors + origin;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Reshape the input vector matrix into its original dimensions.  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_position_vectors                                = reshape(partition_position_vectors, [vector_dimensions indice_dimension_magnitudes]);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%