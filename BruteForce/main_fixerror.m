%% Initialization
clear; close all; clc;

%% (1-) Call the Kinematic Coefficients
% Calculated under KinematicConnections.m

% Loading the Kinematic Coefficients
load( 'KinematicCoefficients.mat' );

% Loading the Robot 
robot = iiwa14( 'high' );
robot.init( );

%% (1A) 3D Case, Position, Hybrid Jacobian, w and w/o Kin. Coefficients, The Wrong Code!

% Get the Symbolic Form of the Hybrid Jacobian
q_syms = sym( 'q', [ 7, 1 ] );
JH = robot.getHybridJacobian( q_syms );

J_sym = JH( 1:3, : );
J_T_sym = J_sym.';

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
JH_val = double( subs( J_sym, q_syms, q_arr ) );

% With Brute-Force Tensor Notation
% a, b are iterating across the joint angles, 1<=a,b<=7
for a = 1 : robot.nq
    for b = 1 : robot.nq
       
        % The Hybrid Jacobian Terms
        % For Translation Only, 1<=i,j<=3
        for i = 1 : 3
            for j = 1 : 3

                % The Derivative Part of J^{i}_a/eta(b)
                diffJ = double( subs( diff( J_sym( i, a ), q_syms( b ) ), q_syms, q_arr ) );

                % Summation
                Kq1( a, b ) = Kq1( a, b ) + JH_val( i, a )*JH_val( j, b )*Kx( i, j ) + diffJ*F_arr( i );

            end
        end
    end
end

Kq1 = double( Kq1 );

% With Johannes Method
% [REF] main_StiffTen_HybridJacobian_3D
% Evaluate and store in matrix

dJ_T_dq = zeros( robot.nq, 3, robot.nq );
dJ_T_dq( :, :, 1 ) = double( subs(  diff( J_T_sym, q_syms( 1 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 2 ) = double( subs(  diff( J_T_sym, q_syms( 2 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 3 ) = double( subs(  diff( J_T_sym, q_syms( 3 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 4 ) = double( subs(  diff( J_T_sym, q_syms( 4 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 5 ) = double( subs(  diff( J_T_sym, q_syms( 5 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 6 ) = double( subs(  diff( J_T_sym, q_syms( 6 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 7 ) = double( subs(  diff( J_T_sym, q_syms( 7 ) ), q_syms, q_arr ) );

K_kin = [ dJ_T_dq( :, :, 1 )*F_arr, dJ_T_dq( :, :, 2 )*F_arr, dJ_T_dq( :, :, 3 )*F_arr, dJ_T_dq( :, :, 4 )*F_arr, ...
          dJ_T_dq( :, :, 5 )*F_arr, dJ_T_dq( :, :, 6 )*F_arr, dJ_T_dq( :, :, 7 )*F_arr ];

J = robot.getHybridJacobian( q_arr );
J = J( 1:3, : );

Kq2 = J' * Kx * J + K_kin;

% The Values for Kq1 and Kq2 are different!!

%% (1B) 3D Case, Position, Hybrid Jacobian, w and w/o Kin. Coefficients, The Correct Code!

% Get the Symbolic Form of the Hybrid Jacobian
q_syms = sym( 'q', [ 7, 1 ] );
JH = robot.getHybridJacobian( q_syms );

J_sym = JH( 1:3, : );
J_T_sym = J_sym.';

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
JH_val = double( subs( J_sym, q_syms, q_arr ) );

% With Brute-Force Tensor Notation
% a, b are iterating across the joint angles, 1<=a,b<=7
for a = 1 : robot.nq
    for b = 1 : robot.nq
       
        % The Hybrid Jacobian Terms
        % For Translation Only, 1<=i,j<=3
        for i = 1 : 3
            for j = 1 : 3

                % Summation
                Kq1( a, b ) = Kq1( a, b ) + JH_val( i, a )*JH_val( j, b )*Kx( i, j );

            end
        end

        % The CODE is FIXED!!
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
% [REF] main_StiffTen_HybridJacobian_3D
% Evaluate and store in matrix

dJ_T_dq = zeros( robot.nq, 3, robot.nq );
dJ_T_dq( :, :, 1 ) = double( subs(  diff( J_T_sym, q_syms( 1 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 2 ) = double( subs(  diff( J_T_sym, q_syms( 2 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 3 ) = double( subs(  diff( J_T_sym, q_syms( 3 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 4 ) = double( subs(  diff( J_T_sym, q_syms( 4 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 5 ) = double( subs(  diff( J_T_sym, q_syms( 5 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 6 ) = double( subs(  diff( J_T_sym, q_syms( 6 ) ), q_syms, q_arr ) );
dJ_T_dq( :, :, 7 ) = double( subs(  diff( J_T_sym, q_syms( 7 ) ), q_syms, q_arr ) );

K_kin = [ dJ_T_dq( :, :, 1 )*F_arr, dJ_T_dq( :, :, 2 )*F_arr, dJ_T_dq( :, :, 3 )*F_arr, dJ_T_dq( :, :, 4 )*F_arr, ...
          dJ_T_dq( :, :, 5 )*F_arr, dJ_T_dq( :, :, 6 )*F_arr, dJ_T_dq( :, :, 7 )*F_arr ];

J = robot.getHybridJacobian( q_arr );
J = J( 1:3, : );

Kq2 = J' * Kx * J + K_kin;

% Kq1 and Kq2 are equivalent!
