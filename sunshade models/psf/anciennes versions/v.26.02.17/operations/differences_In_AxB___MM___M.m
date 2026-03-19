
function    [differences_in_AxB, norms_in_AxB] = ...
			........................................................................................................................................................
            differences_In_AxB___MM___M (set_A, set_B, basis_A, basis_B, origin_A, origin_B)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Assume the first dimension of the sets' matrix representations stores vector information, and the remaining serve as indices.                          %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	List the index dimensions present in each set.   
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


index_dimensions_A = [2:ndims(set_A)];
index_dimensions_B = [2:ndims(set_B)];


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Find the dimensions of the vectors contained within each set (these must be equal).      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


vector_dimensions  = size(set_A, 1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Define the magnitude of the output's first dimensions to correspond to those of set A's indices, and the remaining to those of set B.       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


output_dimensions_differences  = [vector_dimensions size(set_A, index_dimensions_A) size(set_B, index_dimensions_B)];
output_dimensions_norms        = [size(set_A, index_dimensions_A) size(set_B, index_dimensions_B)];


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Calculate the differences between all vector pairs in {AxB}, accouunting for any differences in the origins/bases used in their matrix representations. %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Reshape the matrix representations of the input sets into 2D arrays, such that vectors are expressed along the first dimension.   
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


set_A                   = reshape(set_A, vector_dimensions, []);
set_B                   = reshape(set_B, vector_dimensions, []);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Transform the coordinates of set A and set B vectors (assuming their bases are represented with respect to the same coordinate system).    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


set_A                   = basis_A * set_A;
set_B                   = basis_B * set_B;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Permute the dimensions of set B's representation so that its vectors are listed along its 3rd dimension. Then replicate this 2D array along along the 2nd 
%   dimension to mimick the dimensions of A (so that distances may be taken between each pair in AxB).        
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


set_B                   = permute(set_B,[1 3 2]);
set_B                   = repmat(set_B, [1 length(set_A(1,:)) 1]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Find the differences between each vector pair in AxB, accounting for differences in the origins of the coordinate systems of their matrix representations:
%
%		- Assume that origin_A and origin_B are represented with respect to the same coordinate system.
%		- Last two arguments: take the L2 Norm along the 1st dimension.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%



% PROBLEMS: 

% 	-MULTIPLICATION BY RADII PRIOR TO INPUT NEEDS TO BE CHANGED TO AFTER. 
%	-DOT PRODUCT (A-B)xA HAS TOO MANY DIMENSIONS: NEEDS TO BE DONE BY SHEETS ALONG DIMENSIONS OF A
%	-REORGANIZE CODE HERE AND WITH DOT PRODUCTS FUNCTION TO MAKE MORE EFFICIENT
%	-MORE EFFICIENT NAMING
%	-CHANGE DISTANCES TO NORMS THROUGHOUT AND IN ANALYSIS FUNCTION

differences_in_AxB      = set_A - ((origin_B - origin_A) + set_B);
norms_in_AxB            = vecnorm(differences_in_AxB, 2, 1);
differences_in_AxB      = differences_in_AxB ./ norms_in_AxB;
differences_in_AxB      = reshape(differences_in_AxB, vector_dimensions, []);
differences_in_AxB      = set_A' * differences_in_AxB; 


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Reshape the output matrix so that its first dimensions correspond to the indices of set A and the last to those of set B.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


differences_in_AxB        = reshape(differences_in_AxB, output_dimensions_norms);
norms_in_AxB        = reshape(norms_in_AxB, output_dimensions_norms);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
