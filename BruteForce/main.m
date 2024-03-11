%% (1-) Initialization
clear; close all; clc;

% Loading the Kinematic Coefficients
load( 'KinematicCoefficients.mat' );

% Loading the Robot 
robot = iiwa14( 'high' );
robot.init( );

% Get the Symbolic Form of the Hybrid Jacobian
q_syms = sym( 'q', [ 7, 1 ] );

% The Text 
J_type = [ "Hybrid", "Spatial", "Body" ];

% The J Matrix
J_mat = { robot.getHybridJacobian(  q_syms ), ...
          robot.getSpatialJacobian( q_syms ), ...
          robot.getBodyJacobian( q_syms, 7 ) };

%% (1A) 3D Case, Translation

% Choosing the Type
idx = 3;

fprintf( "\nAnalyzing %s Jacobian, 3D Translation\n", J_type( idx ) );

% Using the Hybrid Jacobian for the Symmetric 
J_sym = J_mat{ idx };
J_T_sym = J_sym.';

% Initializing the joint-stiffness matrix, type 1
Kq1 = sym( zeros( robot.nq, robot.nq ) );

% Initializing the Force, Kx matrices, and Joint array
Kx = diag( [ 402., 311., 167., 0., 0., 0.  ] );

% Check also for non-diagonal matrix
tmp_mat = rand( 3, 3 );
% Kx = tmp_mat * diag( [ 402., 311., 167. ] ) * tmp_mat.';

F_arr = [ 1.5, 2.2, 1.4, 0, 0, 0 ]';
q_arr = robot.q_init;

% The Hybrid Jacobian Value
J_val = double( subs( J_sym, q_syms, q_arr ) );

% With Brute-Force Tensor Notation
% a, b are iterating across the joint angles, 1<=a,b<=7
for a = 1 : robot.nq
    for b = 1 : robot.nq
       
        % The Hybrid Jacobian Terms
        % For Translation Only, 1<=i,j<=3
        for i = 1 : 3
            for j = 1 : 3
                
                tmp = 0;
                if idx == 3 % If Body
                    for k = 1 : 6
                        tmp = tmp + F_arr( k )*Gamma(i,j,k);
                    end
                end

                % Summation
                Kq1( a, b ) = Kq1( a, b ) + J_val( i, a )*J_val( j, b )*( Kx( i, j ) + tmp );

            end
        end

        % The CODE is FIXED!!
        % Must ADD Separately
        tmp = 0;
        for k = 1:3
            diffJ = double( subs( diff( J_sym( k, a ), q_syms( b ) ), q_syms, q_arr ) );
            tmp = tmp + diffJ * F_arr( k );
        end     

        Kq1( a, b ) = Kq1( a, b ) + tmp;

    end
end

Kq1 = double( Kq1 );

% With Johannes Method
% Saving the Jacobian's, nqx6 arrays
dJ_T_dq = zeros( robot.nq, 6, robot.nq );
dJ_T_dq( :, :, 1 ) = double( subs(  diff( J_T_sym, q_syms( 1 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 2 ) = double( subs(  diff( J_T_sym, q_syms( 2 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 3 ) = double( subs(  diff( J_T_sym, q_syms( 3 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 4 ) = double( subs(  diff( J_T_sym, q_syms( 4 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 5 ) = double( subs(  diff( J_T_sym, q_syms( 5 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 6 ) = double( subs(  diff( J_T_sym, q_syms( 6 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 7 ) = double( subs(  diff( J_T_sym, q_syms( 7 ) ), q_syms, q_arr ) );

K_CS = [ 0, 0, 0, 0, -F_arr(3), F_arr(2); ...
         0, 0, 0, F_arr(3), 0, -F_arr(1); ...
         0, 0, 0, -F_arr(2), F_arr(1), 0; ...
         0, 0, 0, 0, 0, 0; ...
         0, 0, 0, 0, 0, 0; ...
         0, 0, 0, 0, 0, 0 ];

K_kin = [ dJ_T_dq( :, :, 1 )*F_arr, dJ_T_dq( :, :, 2 )*F_arr, dJ_T_dq( :, :, 3 )*F_arr, dJ_T_dq( :, :, 4 )*F_arr, ...
          dJ_T_dq( :, :, 5 )*F_arr, dJ_T_dq( :, :, 6 )*F_arr, dJ_T_dq( :, :, 7 )*F_arr ];

if idx == 1 || idx == 2
    Kq2 = J_val' * Kx * J_val + K_kin;   % Without CS

    fprintf( "Error between Brute Force, without CS: %.15f\n", max( max( Kq1 - Kq2  )  ) )
    fprintf( "Symmetry with Brute Force %.15f\n", max( max( Kq2 - Kq2' ) ) );
    fprintf( "Symmetry without CS %.15f\n", max( max( Kq2 - Kq2' ) ) );
    fprintf( "Determinant, Brute Force %.15f\n", det( Kq1 ) );
    fprintf( "Determinant, without CS  %.15f\n", det( Kq2 ) );

elseif idx == 3
    Kq2 = J_val' * ( Kx        ) * J_val + K_kin;   % Without CS    
    Kq3 = J_val' * ( Kx + K_CS ) * J_val + K_kin;   % With    CS

    fprintf( "Error between Brute Force, without CS: %.15f\n", max( max( Kq1 - Kq2  )  ) )
    fprintf( "Error between Brute Force, with    CS: %.15f\n", max( max( Kq1 - Kq3  )  ) )
    fprintf( "Symmetry with Brute Force %.15f\n", max( max( Kq2 - Kq2' ) ) );
    fprintf( "Symmetry without CS %.15f\n", max( max( Kq2 - Kq2' ) ) );
    fprintf( "Symmetry with    CS %.15f\n", max( max( Kq3 - Kq3' ) ) );
    fprintf( "Determinant, Brute Force %.15f\n", det( Kq1 ) );
    fprintf( "Determinant, without CS  %.15f\n", det( Kq2 ) );
    fprintf( "Determinant, with    CS  %.15f\n", det( Kq3 ) );
    
end

%% (1B) 3D Case, Rotation

idx = 2;
fprintf( "\nAnalyzing %s Jacobian, 3D Rotation \n", J_type( idx ) );

% Using the Hybrid Jacobian for the Symmetric 
J_sym = J_mat{ idx };
J_T_sym = J_sym.';

% Initializing the joint-stiffness matrix, type 1
Kq1 = sym( zeros( robot.nq, robot.nq ) );

Kx = diag( [ 0, 0, 0, 40., 31., 16. ] );
F_arr = [ 0, 0, 0, 1.5, 2.2, 1.4 ]';

% Values Randomized vs. Johannes' Code 
q_arr = robot.q_init;

% The Hybrid Jacobian Value
J_val = double( subs( J_sym, q_syms, q_arr ) );

% With Brute-Force Tensor Notation
% a, b are iterating across the joint angles, 1<=a,b<=7
for a = 1 : robot.nq
    for b = 1 : robot.nq
       
        % The Hybrid Jacobian Terms
        % For Rotational Parts

        for i = 4 : 6
            for j = 4 : 6
                tmp = 0;
                if idx == 3 % If Body
                    for k = 4 : 6
                        tmp = tmp + F_arr( k )*Gamma(i,j,k);
                    end
                end

                % Summation
                Kq1( a, b ) = Kq1( a, b ) + J_val( i, a )*J_val( j, b )*( Kx( i, j ) + tmp );

            end
        end

        % The CODE is FIXED!!
        % Must ADD Separately
        tmp = 0;
        for k = 4:6
            diffJ = double( subs( diff( J_sym( k, a ), q_syms( b ) ), q_syms, q_arr ) );
            tmp = tmp + diffJ * F_arr( k );
        end     
        Kq1( a, b ) = Kq1( a, b ) + tmp;

    end
end

Kq1 = double( Kq1 );

% With Johannes Method
% Evaluate and store in matrix
dJ_T_dq = zeros( robot.nq, 6, robot.nq );
dJ_T_dq( :, :, 1 ) = double( subs( diff( J_T_sym, q_syms( 1 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 2 ) = double( subs( diff( J_T_sym, q_syms( 2 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 3 ) = double( subs( diff( J_T_sym, q_syms( 3 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 4 ) = double( subs( diff( J_T_sym, q_syms( 4 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 5 ) = double( subs( diff( J_T_sym, q_syms( 5 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 6 ) = double( subs( diff( J_T_sym, q_syms( 6 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 7 ) = double( subs( diff( J_T_sym, q_syms( 7 ) ), q_syms, q_arr ) );

% Chen 2004!!!!!!!!!!!!!!!!!!!!
K_CS = [ 0, 0, 0, 0, 0, 0; ...
         0, 0, 0, 0, 0, 0; ...
         0, 0, 0, 0, 0, 0; ...
         0, 0, 0, 0, -F_arr(6)/2, F_arr(5)/2; ...
         0, 0, 0, F_arr(6)/2, 0, -F_arr(4)/2; ...
         0, 0, 0, -F_arr(5)/2, F_arr(4)/2, 0 ];

K_kin = [ dJ_T_dq( :, :, 1 ) * F_arr, dJ_T_dq( :, :, 2 ) * F_arr, dJ_T_dq( :, :, 3 ) * F_arr, dJ_T_dq( :, :, 4 ) * F_arr, ...
          dJ_T_dq( :, :, 5 ) * F_arr, dJ_T_dq( :, :, 6 ) * F_arr, dJ_T_dq( :, :, 7 ) * F_arr ];

if idx == 1 || idx == 2
    Kq2 = J_val' * Kx * J_val + K_kin;   % Without CS

    fprintf( "Error between Brute Force, without CS: %.15f\n", max( max( Kq1 - Kq2  )  ) )
    fprintf( "Symmetry with Brute Force %.15f\n", max( max( Kq2 - Kq2' ) ) );
    fprintf( "Symmetry without CS %.15f\n", max( max( Kq2 - Kq2' ) ) );
    fprintf( "Determinant, Brute Force %.15f\n", det( Kq1 ) );
    fprintf( "Determinant, without CS  %.15f\n", det( Kq2 ) );

elseif idx == 3
    Kq2 = J_val' * ( Kx        ) * J_val + K_kin;   % Without CS    
    Kq3 = J_val' * ( Kx + K_CS ) * J_val + K_kin;   % With    CS

    fprintf( "Error between Brute Force, without CS: %.15f\n", max( max( Kq1 - Kq2  )  ) )
    fprintf( "Error between Brute Force, with    CS: %.15f\n", max( max( Kq1 - Kq3  )  ) )
    fprintf( "Symmetry with Brute Force %.15f\n", max( max( Kq2 - Kq2' ) ) );
    fprintf( "Symmetry without CS %.15f\n", max( max( Kq2 - Kq2' ) ) );
    fprintf( "Symmetry with    CS %.15f\n", max( max( Kq3 - Kq3' ) ) );
    fprintf( "Determinant, Brute Force %.15f\n", det( Kq1 ) );
    fprintf( "Determinant, without CS  %.15f\n", det( Kq2 ) );
    fprintf( "Determinant, with    CS  %.15f\n", det( Kq3 ) );
    
end
%% (1C) 6D Case Translation and Rotation

idx = 3;
fprintf( "\nAnalyzing %s Jacobian, 6D Translation + Rotation \n", J_type( idx ) );

% Using the Hybrid Jacobian for the Symmetric 
J_sym = J_mat{ idx };
J_T_sym = J_sym.';

% Initializing the joint-stiffness matrix, type 1
Kq1 = sym( zeros( robot.nq, robot.nq ) );

Kx = diag( [ 402., 311., 167., 42., 31., 17. ] );
F_arr = [ 1.5, 2.2, 1.4, 0.7, 1.2, 0.95 ]';

% Values Randomized vs. Johannes' Code 
q_arr = robot.q_init;

% The Hybrid Jacobian Value
J_val = double( subs( J_sym, q_syms, q_arr ) );

% With Brute-Force Tensor Notation
% a, b are iterating across the joint angles, 1<=a,b<=7
for a = 1 : robot.nq
    for b = 1 : robot.nq
       
        % The Hybrid Jacobian Terms
        % For Rotational Parts
        for i = 1 : 6
            for j = 1 : 6
                
                tmp = 0;

                if idx == 3 % If Body
                    for k = 1 : 6
                        tmp = tmp + F_arr( k )*Gamma(i,j,k);
                    end
                end

                % Summation
                Kq1( a, b ) = Kq1( a, b ) + J_val( i, a )*J_val( j, b )*( Kx( i, j ) + tmp );

            end
        end

        % The CODE is FIXED!!
        % Must ADD Separately
        tmp = 0;
        for k = 1:6
            diffJ = double( subs( diff( J_sym( k, a ), q_syms( b ) ), q_syms, q_arr ) );
            tmp = tmp + diffJ * F_arr( k );
        end     
        Kq1( a, b ) = Kq1( a, b ) + tmp;


    end
end

Kq1 = double( Kq1 );

% With Johannes MethodW
% Evaluate and store in matrix
dJ_T_dq = zeros( robot.nq, 6, robot.nq );
dJ_T_dq( :, :, 1 ) = double( subs( diff( J_T_sym, q_syms( 1 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 2 ) = double( subs( diff( J_T_sym, q_syms( 2 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 3 ) = double( subs( diff( J_T_sym, q_syms( 3 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 4 ) = double( subs( diff( J_T_sym, q_syms( 4 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 5 ) = double( subs( diff( J_T_sym, q_syms( 5 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 6 ) = double( subs( diff( J_T_sym, q_syms( 6 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 7 ) = double( subs( diff( J_T_sym, q_syms( 7 ) ), q_syms, q_arr ) );


% For 6D
K_CS1 = [ 0, 0, 0, 0, -F_arr(3), F_arr(2); ...
        0, 0, 0, F_arr(3), 0, -F_arr(1); ...
        0, 0, 0, -F_arr(2), F_arr(1), 0; ...
        0, 0, 0, 0, -F_arr(6)/2, F_arr(5)/2; ...
        0, 0, 0, F_arr(6)/2, 0, -F_arr(4)/2; ...
        0, 0, 0, -F_arr(5)/2, F_arr(4)/2, 0 ];

% Chen 2004 !!!!!!!!
% This is wrong with the full tensor notation
K_CS2 = [ 0, 0, 0, 0, 0, 0; ...
          0, 0, 0, 0, 0, 0; ...
          0, 0, 0, 0, 0, 0; ...
          0, 0, 0, 0, -F_arr(6)/2, F_arr(5)/2; ...
          0, 0, 0,  F_arr(6)/2, 0, -F_arr(4)/2; ...
          0, 0, 0, -F_arr(5)/2, F_arr(4)/2, 0 ];

K_kin = [ dJ_T_dq( :, :, 1 ) * F_arr, dJ_T_dq( :, :, 2 ) * F_arr, dJ_T_dq( :, :, 3 ) * F_arr, dJ_T_dq( :, :, 4 ) * F_arr, ...
          dJ_T_dq( :, :, 5 ) * F_arr, dJ_T_dq( :, :, 6 ) * F_arr, dJ_T_dq( :, :, 7 ) * F_arr ];

if idx == 3 % If Body
    Kq2 = J_val' * Kx * J_val + J_val' * K_CS1 * J_val + K_kin;
    Kq3 = J_val' * Kx * J_val + J_val' * K_CS2 * J_val + K_kin;
    Kq4 = J_val' * Kx * J_val + K_kin;    
else
    Kq2 = J_val' * Kx * J_val + K_kin;
    Kq3 = J_val' * Kx * J_val + K_kin;
end

% Kq1 and Kq2 are equivalent!
fprintf( "Error between two methods with    off-diag: %.15f\n", max( max( Kq1 - Kq2  )  ) )
fprintf( "Error between two methods without off-diag: %.15f\n", max( max( Kq1 - Kq3  )  ) )

% Symmetry
fprintf( "Symmetry with    off-diag %.15f\n", max( max( Kq2 - Kq2' ) ) );
fprintf( "Symmetry without off-diag %.15f\n", max( max( Kq3 - Kq3' ) ) );
fprintf( "Determinant %.15f\n", det( Kq1 ) );

if idx == 3
    fprintf( "Symmetry only Kin Stiff %.15f\n", max( max( Kq4 - Kq4' ) ) );
end
