
function    delta_deg = ...
            ........................................................................................................................................................
            spencer_declination (day_of_year, days_in_year)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Solar declination angle computed via the Spencer (1971) Fourier series.                                                                                  %&%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Reference:
%       Spencer, J.W. (1971). "Fourier series representation of the position
%       of the sun." Search, 2(5), 172.
%
%   The formula requires no external data and is valid for forward projection
%   to any future date, making it appropriate for climate modelling scenarios.
%
%   Inputs:
%       day_of_year   - scalar or array; day number within the year (1 = Jan 1)
%       days_in_year  - (optional) 365 for a regular year, 366 for a leap year.
%                       Defaults to 365 if omitted.
%
%   Output:
%       delta_deg     - solar declination angle in degrees (same size as input).
%                       Positive values indicate the sub-stellar point is north
%                       of the geographic equator (northern summer).
%
%   Accuracy: within ~0.035 deg of the exact value across the full year.
%
%   Usage in coordinate correction:
%       The sunshade irradiance model computes shading in the heliocentric
%       frame, where phi = 0 always coincides with the sub-stellar point.
%       To map results to geographic latitude, the output array is shifted
%       by round(delta_deg / lat_step_deg) rows before writing to the NC file.
%       See analysis_General___0___X.m for implementation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Default year length.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


if nargin < 2
    days_in_year = 365;
end


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Day-angle B in radians: fraction of the full annual cycle elapsed by day d.
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


B = (2 * pi / days_in_year) .* (day_of_year - 1);


%------------------------------------------------------------------------------------------------------------------------------------------------------------------%
%   Spencer (1971) Fourier series for declination (result in radians, converted to degrees).
%------------------------------------------------------------------------------------------------------------------------------------------------------------------%


delta_deg = (180 / pi) .* ( ...
      0.006918 ...
    - 0.399912 .* cos(     B) ...
    + 0.070257 .* sin(     B) ...
    - 0.006758 .* cos( 2 * B) ...
    + 0.000907 .* sin( 2 * B) ...
    - 0.002697 .* cos( 3 * B) ...
    + 0.00148  .* sin( 3 * B) );


%%&%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
