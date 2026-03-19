
function    [dot_products_in_AxB] = ...
			........................................................................................................................................................
            dot_products_In_AxB___4M___M (set_A, set_B, basis_A, basis_B)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Calculate the dot product of each vector pair in AxB and output the results with respect to the original vector indices.                                %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Record the magnitude of the indice dimensions of the input matrix representations (all but the 1st dimension, which is assumed to store vector information),
%   reshape the input set of vectors into a 2D matrix, and perform changes of coordinates to the same basis (via input basis representations).       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[set_A, indice_dimension_magnitudes_A] = represent_Vectors_As_2D_Matrix___2M___MR (set_A, basis_A);
[set_B, indice_dimension_magnitudes_B] = represent_Vectors_As_2D_Matrix___2M___MR (set_B, basis_B);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the dot products of all vector pairs in {AxB} via matrix multiplication.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_AxB                    = set_A' * set_B;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Reshape the output matrix so that its first dimensions correspond to the indices of set A and the last to those of set B.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_AxB                    = reshape(dot_products_in_AxB, [indice_dimension_magnitudes_A indice_dimension_magnitudes_B]);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
