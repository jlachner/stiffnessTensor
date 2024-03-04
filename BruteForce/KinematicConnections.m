%% Initialization
clear; close all; clc;

%% Calculating the Connection Coefficients

% Calculating the Structure Constants
e1 = [ zeros( 3, 3 ), [1;0;0]; zeros(1,4 )];
e2 = [ zeros( 3, 3 ), [0;1;0]; zeros(1,4 )];
e3 = [ zeros( 3, 3 ), [0;0;1]; zeros(1,4 )];

e4 = [ 0,  0,  0, 0; 
       0,  0, -1, 0;
       0,  1,  0, 0;
       zeros( 1, 4 )];

e5 = [ 0,  0, 1, 0; 
       0,  0, 0, 0;
      -1,  0, 0, 0;
       zeros( 1, 4 )];

e6 = [ 0, -1, 0, 0; 
       1,  0, 0, 0;
       0,  0, 0, 0;
       zeros( 1, 4 )];

e_arr = zeros( 4, 4, 6 );
e_arr( :, :, 1 ) = e1;
e_arr( :, :, 2 ) = e2;
e_arr( :, :, 3 ) = e3;
e_arr( :, :, 4 ) = e4;
e_arr( :, :, 5 ) = e5;
e_arr( :, :, 6 ) = e6;

v1 = se3_to_R6( e1 );
v2 = se3_to_R6( e2 );
v3 = se3_to_R6( e3 );
v4 = se3_to_R6( e4 );
v5 = se3_to_R6( e5 );
v6 = se3_to_R6( e6 );

v_arr = zeros( 6, 6 );
v_arr( :, 1 ) = v1;
v_arr( :, 2 ) = v2;
v_arr( :, 3 ) = v3;
v_arr( :, 4 ) = v4;
v_arr( :, 5 ) = v5;
v_arr( :, 6 ) = v6;

% Construction coefficients
Cijk = zeros( 6, 6, 6 );

% Calculation
for i = 1 : 6
    for j = 1 : 6
        tmp1 = e_arr( :, :, i )*e_arr( :, :, j )-e_arr( :, :, j )*e_arr( :, :, i );
        tmp2 = se3_to_R6( tmp1 );        
        for k = 1: 6
            val = sum( tmp2 .* v_arr( :, k )' );

            if val
                fprintf( "C_{%d, %d}^{%d} = %.2f\n", i, j, k, val )
            end
            Cijk( i, j, k ) = val;
        end
    end
end

% Also calculate the a values.
% Manual Calculation.
aijk = zeros( 6, 6, 6);
aijk( 4, 2, 3 ) = +1;
aijk( 4, 3, 2 ) = -1;
aijk( 5, 1, 3 ) = -1;
aijk( 6, 1, 2 ) = +1;
aijk( 5, 3, 1 ) = +1;
aijk( 6, 2, 1 ) = -1;


% Calculating Gammas
Gamma = zeros( 6, 6, 6 );
for i = 1:6
    for j = 1:6
        for k = 1:6
            if aijk( i, j , k ) || Cijk( i, j, k ) 
                    Gamma( j, i, k ) = 0.5 * ( +Cijk( i, j, k  ) + aijk( i, j , k ) );
                    Gamma( i, j, k ) = 0.5 * ( -Cijk( i, j, k  ) + aijk( i, j , k ) );
            end
        end
    end
end


for i = 1 : 6
    for j = 1 : 6   
        for k = 1: 6
            if Gamma( i, j, k )
                fprintf( "Gamma_{%d, %d}^{%d} = %.2f\n", i, j, k, Gamma( i, j, k ) )
            end
        end
    end
end

% Deriving the force
syms F1 F2 F3 F4 F5 F6 
F_arr = [ F1; F2; F3; F4; F5; F6 ];

Kasym = sym( zeros( 6, 6 ) );

for i = 1 : 6
    for j = 1 : 6
        tmp = 0;
        for k = 1 : 6
            tmp = tmp + Gamma( i,j, k ) * F_arr( k );
        end
        Kasym( i, j ) = tmp;
    end
end

% Save the Gamma Coefficients
save( 'KinematicCoefficients.mat', 'Gamma')