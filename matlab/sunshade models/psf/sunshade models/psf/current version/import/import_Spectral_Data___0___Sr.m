
function    [spectral_data] = ...
			........................................................................................................................................................
            import_Spectral_Data___0___Sr


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Import the spectral data and perform .                                  					%&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the Excel file location of the spectral data and its corresponding locations within the file.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


[excel_file]                                                  = location_Star_Spectral_Data___E___Sr;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Import the spectral data.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


spectral_data_import                                          = readmatrix(excel_file.location, 'Range', excel_file.range.all_data);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Assign structure variables to spectral data elements.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


spectral_data.bands                                           = spectral_data_import(:,excel_file.columns.spectral_bands);
spectral_data.irradiances                                     = spectral_data_import(:,excel_file.columns.spectral_irradiances);
spectral_data.coefficients                                    = spectral_data_import(:,excel_file.columns.spectral_coefficients);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Convert spectral irradiance and coefficient values imported as NaN (i.e. no data entered in Excel file) to zero.      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


spectral_data.irradiances(isnan(spectral_data.irradiances))   = 0;
spectral_data.coefficients(isnan(spectral_data.coefficients)) = 0;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Normalize the irradiances amongst the spectral bands by attributing all non-assigned irradiance to the default band.   
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


spectral_data.irradiances(1,1)                                = 1-sum(spectral_data.irradiances(2:end,1));


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Count the number of spectral bands.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


spectral_data.number_of_bands                                 = size(spectral_data.bands,1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Calculate the normalizing constant for each spectral band (each coefficient's contribution is calculated one at a time, then divided by a constant).      
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


spectral_data.normalizing_coefficients                        = 8*spectral_data.coefficients(:,4);
spectral_data.normalizing_coefficients                        = spectral_data.normalizing_coefficients + 3*pi*spectral_data.coefficients(:,3);
spectral_data.normalizing_coefficients                        = spectral_data.normalizing_coefficients + 12*spectral_data.coefficients(:,2);
spectral_data.normalizing_coefficients                        = spectral_data.normalizing_coefficients + 6*pi*spectral_data.coefficients(:,1);
spectral_data.normalizing_coefficients                        = spectral_data.normalizing_coefficients / 6;


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%	Normalize the spectral coefficients.     
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


spectral_data.coefficients                                    = spectral_data.coefficients ./ spectral_data.normalizing_coefficients;


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%