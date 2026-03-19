
function    [output_matrix] = ...
			........................................................................................................................................................
            verification_Distances_In_AxB___MM___M (input_matrix_A, input_matrix_B, distance_between_origins)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Calculate the distances between all pairs of vectors in {AxB} (expressed along the 1st dimension of the respective input matrices).                     %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Find the second and third dimensions of the input matrices.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


input_matrix_A_dim_2                   = length(input_matrix_A(1,:,1));
input_matrix_A_dim_3                   = length(input_matrix_A(1,1,:));

input_matrix_B_dim_2                   = length(input_matrix_B(1,:,1));
input_matrix_B_dim_3                   = length(input_matrix_B(1,1,:));


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Express the distance between the origins of the 2 vector sets as a vector along the x-axis, from the vector set with positive x-values to that with negative.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if max(input_matrix_A(:,:,1)) >= 0

	vector_between_origins             = [distance_between_origins ; 0 ; 0]

else

	vector_between_origins             = [-distance_between_origins ; 0 ; 0]

end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Assign the first two dimensions of the output matrix as those of input matrix A, and the remaining two as those of input matrix B. 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


output_matrix                          = zeros(input_matrix_A_dim_2, input_matrix_A_dim_3, input_matrix_B_dim_2, input_matrix_B_dim_3);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the distance between each possible pair of vectors in {AxB} (expressed along the 3rd dimension of the respective input matrices), noting:    
%
%		- the x-axis of the vectors defined in set A is oriented opposite to the x-asix of those defined in set B
%		- 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------% 


for i = 1:input_matrix_A_dim_2

	for j = 1:input_matrix_A_dim_3

		for n = 1:input_matrix_B_dim_2

			for m = 1:input_matrix_B_dim_3

				output_matrix(i,j,n,m) = norm(input_matrix_A(:,i,j) - vector_between_origins - input_matrix_B(:,n,m));

			end

		end

	end

end
				

%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%