
function    [output_matrix] = ...
			........................................................................................................................................................
            verification_Dot_Products_In_AxB___MM___M (input_matrix_A, input_matrix_B)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Calculate the dot products of all possible pairs of vectors in {AxB} (expressed along the 1st dimension of the respective input matrices).              %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Find the second and third dimensions of the input matrices.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


input_matrix_A_dim_2                   = length(input_matrix_A(1,:,1));
input_matrix_A_dim_3                   = length(input_matrix_A(1,1,:));

input_matrix_B_dim_2                   = length(input_matrix_B(1,:,1));
input_matrix_B_dim_3                   = length(input_matrix_B(1,1,:));


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Assign the first two dimensions of the output matrix as those of input matrix A, and the remaining two as those of input matrix B. 
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


output_matrix                          = zeros(input_matrix_A_dim_2, input_matrix_A_dim_3, input_matrix_B_dim_2, input_matrix_B_dim_3);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the dot products of each possible pair of vectors in {AxB} (expressed along the 1st dimension of the respective input matrices).     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------% 


for i = 1:input_matrix_A_dim_2

	for j = 1:input_matrix_A_dim_3

		for n = 1:input_matrix_B_dim_2

			for m = 1:input_matrix_B_dim_3

				output_matrix(i,j,n,m) = dot(input_matrix_A(:,i,j), input_matrix_B(:,n,m));

			end

		end

	end

end
				

%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%