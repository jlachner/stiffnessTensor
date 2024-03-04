%% Initialization
clear; close all; clc;

%% (1-) Call the Kinematic Coefficients
% Calculated under KinematicConnections.m

% Loading the Kinematic Coefficients
load( 'KinematicCoefficients.mat' );

% Loading the Robot 
robot = iiwa14( 'high' );
robot.init( );

%% (1A) 3D Case, Position, Hybrid Jacobian, w and w/o Kin. Coefficients

% Get the Symbolic Form of the Hybrid Jacobian
q_syms = sym( 'q', [ 7, 1 ] );
JH = robot.getHybridJacobian( q_syms );

% Initializing the joint-stiffness matrix, type 1
Kq1 = sym( zeros( robot.nq, robot.nq ) );

% Initializing the Force, Kx matrices, and Joint array
% Randomized
% tmp1 = rand( 3, 3 );
% Kx = tmp1' * diag( [ 1.2, 4.3, 5.4 ] ) * tmp1;

Kx = zeros( 3, 3 );
Kx( 1, 1 ) = 402;
Kx( 2, 2 ) = 311;
Kx( 3, 3 ) = 167;

% Values Randomized vs. Johannes' Code 
% F_arr = rand( 3, 1 );
F_arr = [ 1.5, 2.2, 1.4 ]';

% Values Randomized vs. Johannes' Code 
% q_arr = rand( robot.nq,1 );
q_arr = robot.q_init;

% The Hybrid Jacobian Value
JH_val = double( subs( JH, q_syms, q_arr ) );

% With Brute-Force Tensor Notation
% a, b are iterating across the joint angles, 1<=a,b<=7
for a = 1 : robot.nq
    for b = 1 : robot.nq
        
        K_ab = 0;

        % The Hybrid Jacobian Terms
        % For Translation Only, 1<=i,j<=3
        for i = 1 : 3
            for j = 1 : 3
                tmp = 0;

                % Note that the translational term has no connection coefficients
                % So the following term is not needed, strictly speaking
                for k = 1 : 3
                    tmp = tmp + F_arr( k )*Gamma(i,j,k);
                end

                % The Derivative Part of J^{i}_a/eta(b)
                diffJ = double( subs( diff( JH( i, a ), q_syms( b ) ), q_syms, q_arr ) );

                % Summation
                K_ab = K_ab + JH_val( i, a )*JH_val( j, b )*( Kx( i, j ) + tmp ) + diffJ * F_arr( i );

            end
        end

        Kq1( a, b ) = K_ab;
    end
end

Kq1 = double( Kq1 );

% With Johannes Method
% [REF] main_StiffTen_HybridJacobian_3D
J_sym = JH( 1:3, : );
J_T_sym = J_sym.';

% Partial derivative of J_T_sym wrt. q_sym
dJ_T_dq1 = diff( J_T_sym, q_syms( 1 ) );
dJ_T_dq2 = diff( J_T_sym, q_syms( 2 ) );
dJ_T_dq3 = diff( J_T_sym, q_syms( 3 ) );
dJ_T_dq4 = diff( J_T_sym, q_syms( 4 ) );
dJ_T_dq5 = diff( J_T_sym, q_syms( 5 ) );
dJ_T_dq6 = diff( J_T_sym, q_syms( 6 ) );
dJ_T_dq7 = diff( J_T_sym, q_syms( 7 ) );

% Evaluate and store in matrix
dJ_T_dq = zeros( robot.nq, 3, robot.nq );
dJ_T_dq( :, :, 1 ) = double( subs( dJ_T_dq1, q_syms, q_arr ) );
dJ_T_dq( :, :, 2 ) = double( subs( dJ_T_dq2, q_syms, q_arr ) );
dJ_T_dq( :, :, 3 ) = double( subs( dJ_T_dq3, q_syms, q_arr ) );
dJ_T_dq( :, :, 4 ) = double( subs( dJ_T_dq4, q_syms, q_arr ) );
dJ_T_dq( :, :, 5 ) = double( subs( dJ_T_dq5, q_syms, q_arr ) );
dJ_T_dq( :, :, 6 ) = double( subs( dJ_T_dq6, q_syms, q_arr ) );
dJ_T_dq( :, :, 7 ) = double( subs( dJ_T_dq7, q_syms, q_arr ) );

K_kin = [ dJ_T_dq( :, :, 1 ) * F_arr, dJ_T_dq( :, :, 2 ) * F_arr, dJ_T_dq( :, :, 3 ) * F_arr, dJ_T_dq( :, :, 4 ) * F_arr, ...
          dJ_T_dq( :, :, 5 ) * F_arr, dJ_T_dq( :, :, 6 ) * F_arr, dJ_T_dq( :, :, 7 ) * F_arr ];

J = robot.getHybridJacobian( q_arr );
J = J( 1:3, : );

Kq2 = J' * Kx * J + K_kin;
