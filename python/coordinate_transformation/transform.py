from coordinate_transformation.assemble import assemble_from_files

assemble_from_files(
    source_nc_path="/Users/morgangoodwin/Desktop/PSI/Matlab/numerical control/exports/L1-RX-2026-04-16-001.nc",
    horizons_txt_path="/Users/morgangoodwin/Desktop/PSI/Matlab/Planning/Coordinate_Transformation/horizons_Sun_RA-DEC.txt",
    output_nc_path="L1-RX-2026-04-16-001_df_2035_2065.nc",
    show_progress=True,
)
