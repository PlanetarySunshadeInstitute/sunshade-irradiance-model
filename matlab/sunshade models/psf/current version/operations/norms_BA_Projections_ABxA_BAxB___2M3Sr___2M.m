
function    [norms_BA, projections_ABxA, projections_BAxB ] = ...
			........................................................................................................................................................
            norms_BA_Projections_ABxA_BAxB___2M3Sr___2M (sets_of_normal_vectors, sets_of_position_vectors, bases, origins)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Identify the bases, origins, and scaling factors corresponding to each set assuming thesse were each input as 2D structures.                            %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Identify the sets, bases, origins, and scaling factor structure fieldnames and define these as individual variables.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


set_fieldnames             = fieldnames(sets_of_normal_vectors);
set_normals_A              = sets_of_normal_vectors.(set_fieldnames{1});
set_normals_B              = sets_of_normal_vectors.(set_fieldnames{2});


set_fieldnames             = fieldnames(sets_of_position_vectors);
set_positions_A            = sets_of_position_vectors.(set_fieldnames{1});
set_positions_B            = sets_of_position_vectors.(set_fieldnames{2});


bases_fieldnames           = fieldnames(bases);
basis_A                    = bases.(bases_fieldnames{1});
basis_B                    = bases.(bases_fieldnames{2});


origins_fieldnames         = fieldnames(origins);
origin_A                   = origins.(origins_fieldnames{1});
origin_B                   = origins.(origins_fieldnames{2});


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Clear the original sets of vectors for efficiency.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


clear sets_of_normal_vectors


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Assume the first dimension of the sets' matrix representations stores vector information, and the remaining serve as indices.                          %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Record the magnitude of the indice dimensions of the input matrix representations (all but the 1st dimension, which is assumed to store vector information),
%   reshape the input sets of vectors into 2D matrices, and perform changes of coordinates to the same basis (assumed via input basis representations).       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[set_normals_A, indice_dimension_magnitudes_A]   = represent_Vectors_As_2D_Matrix___2M___MR (set_normals_A, basis_A);
[set_normals_B, indice_dimension_magnitudes_B]   = represent_Vectors_As_2D_Matrix___2M___MR (set_normals_B, basis_B);

[set_positions_A, ~]                             = represent_Vectors_As_2D_Matrix___2M___MR (set_positions_A, basis_A);
[set_positions_B, ~]                             = represent_Vectors_As_2D_Matrix___2M___MR (set_positions_B, basis_B);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Permute the dimensions of the 2D array of set B's vectors so that they're now listed along its (once trivial) 3rd dimension. Then replicate the resulting
%   2D (dimensions 1,3) sub-array along the 2nd dimension to correspond to the magnitude of the 2nd dimensions of (number of vectors in) set A.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


set_normals_B                                    = permute(set_normals_B,[1 3 2]);
set_normals_B                                    = repmat(set_normals_B, [1 length(set_normals_A(1,:)) 1]);


set_positions_B                                  = permute(set_positions_B,[1 3 2]);
set_positions_B                                  = repmat(set_positions_B, [1 length(set_positions_A(1,:)) 1]);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Take the differences between all vectors in sets A and B, multiplying by respective scaling factors and accounting for differences in the origins of their 
%	input representations.   
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


set_BA                                           = set_positions_B - ((origin_A - origin_B) + set_positions_A);                                  


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the norms of these differences (last two arguments: take the L2 norm along the 1st dimension), then normalize the difference vectors.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


norms_BA                                         = vecnorm(set_BA, 2, 1);
set_BA                                           = set_BA ./ norms_BA;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Begin taking the dot products in ABxA by performing element-wise multiplication (of each Axb slice of the 3D array of differences in AxB by set A's 2D array,
%	and similarly with set B as it has been reformatted above).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


projections_ABxA                                 = set_BA .* set_normals_A;
projections_BAxB                                 = -set_BA .* set_normals_B;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	For efficiency, clear the explicit information defining the set of differences in AxB and the vectors in sets A and B, which is no longer needed.    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


clear                                              set_BA set_normals_A set_normals_B


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Finish taking the dot products by summing their element-wise products over the first dimension, then composing with the indicator function on the positive
%	reals [whereby only angles of incidence in (-90,90) are counted].  
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


projections_ABxA                                 = sum(projections_ABxA, 1);
projections_ABxA                                 = max(0, projections_ABxA);

projections_BAxB                                 = sum(projections_BAxB, 1);
projections_BAxB                                 = max(0, projections_BAxB);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Reshape the arrays of scalar results so that the first dimensions correspond to the indices of set A and the latter to those of set B.       
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


projections_ABxA                                 = reshape(projections_ABxA, [indice_dimension_magnitudes_A indice_dimension_magnitudes_B]);
projections_BAxB                                 = reshape(projections_BAxB, [indice_dimension_magnitudes_A indice_dimension_magnitudes_B]);
norms_BA                                         = reshape(norms_BA, [indice_dimension_magnitudes_A indice_dimension_magnitudes_B]);


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
