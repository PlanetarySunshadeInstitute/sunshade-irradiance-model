
function    [partition_vectors] = ...
			........................................................................................................................................................
            partition_Scale_Transform_Shift___2MCS___Sr (partition_vectors, change_of_basis, origin, scaling_factor)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Perform a change of basis on the position and normal vectors, relocate the origin of the position vectors, and normalize the normal vectors.            %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Generate the identity matrix corresponding to the dimensions of the input vectors (the magnitude of the first dimension of their input matrix), record this
%	matrix's initial dimensions, then reshape it into a 2D matrix (where vector information is stored along the first dimension).      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


identity_transformation                          = eye(size(partition_vectors,1));
[partition_vectors, indice_dimension_magnitudes] = represent_Vectors_As_2D_Matrix___2M___MR (partition_vectors, identity_transformation);
vector_dimensions                                = size(partition_vectors,1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Scale the partition's vectors by the input scaling factor.       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_vectors = scaling_factor * partition_vectors;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Perform the input change of basis transformation.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_vectors = change_of_basis * partition_vectors;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Shift the origin of the partition's coordinate system to that of the common.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_vectors = partition_vectors + origin;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Reshape the input vector matrix into its original dimensions.  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


partition_vectors = reshape(partition_vectors, [vector_dimensions indice_dimension_magnitudes]);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%