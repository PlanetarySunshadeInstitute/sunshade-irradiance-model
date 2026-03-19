
function    [norms_in_AxB, dot_products_in_ABxA] = ...
			........................................................................................................................................................
            norms_In_AxB_And_Dot_Products_In_ABxA___4M2C2S___2M (set_A, set_B, basis_A, basis_B, origin_A, origin_B, scaling_factor_A, scaling_factor_B)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Assume the first dimension of the sets' matrix representations stores vector information, and the remaining serve as indices.                          %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Record the magnitude of the indice dimensions of the input matrix representations (all but the 1st dimension, which is assumed to store vector information),
%   reshape the input set of vectors into a 2D matrix, and perform changes of coordinates to the same basis (via input basis representations).       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[set_A, indice_dimension_magnitudes_A] = represent_Vectors_As_2D_Matrix___2M___MR (set_A, basis_A);
[set_B, indice_dimension_magnitudes_B] = represent_Vectors_As_2D_Matrix___2M___MR (set_B, basis_B);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Permute the dimensions of the 2D array of set B's vectors so that they're now listed along its (once trivial) 3rd dimension. Then replicate the resulting
%   2D (dimensions 1,3) sub-array along the 2nd dimension to correspond to the magnitude of the 2nd dimensions of (number of vectors in) set A.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


set_B                                  = permute(set_B,[1 3 2]);
set_B                                  = repmat(set_B, [1 length(set_A(1,:)) 1]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Take the differences between all vectors in sets A and B, multiplying by respective scaling factors and accounting for differences in the origins of their 
%	input representations.   
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


differences_in_AxB                     = scaling_factor_B * set_B - ((origin_A - origin_B) + scaling_factor_A * set_A);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	For efficiency, clear the vectors in set B which are no longer needed.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


clear                                    set_B


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the norms of these differences (last two arguments: take the L2 norm along the 1st dimension), then normalize the difference vectors.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


norms_in_AxB                           = vecnorm(differences_in_AxB, 2, 1);
differences_in_AxB                     = differences_in_AxB ./ norms_in_AxB;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Begin taking the dot products in ABxA by performing element-wise multiplication (of each Axb slice of the 3D array of differences in AxB by set A's 2D array.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_ABxA                   = differences_in_AxB .* set_A;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	For efficiency, clear the set of differences in AxB and the vectors in set A, which are no longer needed.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


clear                                    differences_in_AxB set_A


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Finish taking the dot products in ABxA by summing their element-wise products over the first dimension, then composing with the indicator function on the 
%	positive reals (this last is redundant via multiplication by dot products in AxB, but can reduce data size).  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_ABxA                   = sum(dot_products_in_ABxA, 1);
dot_products_in_ABxA                   = max(0, dot_products_in_ABxA);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Reshape the arrays of scalar results so that the first dimensions correspond to the indices of set A and the latter to those of set B.       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


dot_products_in_ABxA                   = reshape(dot_products_in_ABxA, [indice_dimension_magnitudes_A indice_dimension_magnitudes_B]);
norms_in_AxB                           = reshape(norms_in_AxB, [indice_dimension_magnitudes_A indice_dimension_magnitudes_B]);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
