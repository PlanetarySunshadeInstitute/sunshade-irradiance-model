
function    [output_matrix] = ...
			........................................................................................................................................................
            dot_products_of_all_3rd_dimension_combinations___MM___M (input_matrix_i, input_matrix_f)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Calculate the dot products of all possible combinations of the vectors expressed along the 3rd dimension of the input matrices.                         %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Find the first two dimensions of the input matrices.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


input_matrix_i_dim_1                   = length(input_matrix_i(:,1,1));
input_matrix_i_dim_2                   = length(input_matrix_i(1,:,1));

input_matrix_f_dim_1                   = length(input_matrix_f(:,1,1));
input_matrix_f_dim_2                   = length(input_matrix_f(1,:,1));


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Assign the first two dimensions of the output matrix as those of the final input matrix, and the remaining two as those of the initial input matrix. 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


output_matrix                          = zeros(input_matrix_f_dim_1, input_matrix_f_dim_2, input_matrix_i_dim_1, input_matrix_i_dim_2);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the dot product of each possible combination of vectors expressed along the 3rd dimension of the input matrix.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------% 


for i = 1:input_matrix_f_dim_1

	for j = 1:input_matrix_f_dim_2

		for n = 1: input_matrix_i_dim_1

			for m = 1:input_matrix_i_dim_2

				output_matrix(i,j,k,l) = dot(input_matrix_f(i,j,:), input_matrix_i(n,m,:));

			end

		end

	end

end
				

%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%