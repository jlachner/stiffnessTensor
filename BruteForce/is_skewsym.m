function is_check = is_skewsym( mat )
% ===========================================================================
% Descriptions
% ------------
%    Check whether the matrix is asymmetric, i.e., an element of so(3)
% 
% Parameters
% ----------
%   (1)   mat: 3x3 matrix
%   (2) thres: scalar
% 
% Returns
% -------
%   (1) is_check: boolean value of true/false
%
% ===========================================================================

% Should be a 3x3 matrix and positive values
assert( all( size( mat ) == [ 3, 3 ] ) );
% assert( thres > 0 );

% Get the error value with norm
err = norm( mat + mat' );

% Threshold value for evaluating so3
thres = 1e-5;

% Return whether true ro false
is_check = ( err <= 1e-7 );

end