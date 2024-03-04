% [Project]        Robot Simulator - iiwa14
% Authors                       Email
%   [1] Johannes Lachner        jlachner@mit.edu
%   [2] Moses C. Nah            mosesnah@mit.edu
%
%
% The code is heavily commented. A famous quote says:
% "Code is read more often than it is written"
%           - Guido Van Rossum, the Creator of Python

%% Cleaning up + Environment Setup
clear; close all; clc;

% Simulation settings
simTime = 5;        % Total simulation time
t  = 0;             % The current time of simulation   
dt = 0.01;          % Time-step of simulation 

% Set figure size and attach robot to simulation
robot = iiwa14( 'low' );
nq = robot.nq;
robot.init( );

%% Create animation
anim = Animation( 'Dimension', 3, 'xLim', [-0.7,0.7], 'yLim', [-0.7,0.7], 'zLim', [0,1.4] );
% anim.init( );
% anim.attachRobot( robot )  

%% Update kinematics
q = robot.q_init;
robot.updateKinematics( q );
% anim.update( 0 );

%% Kinematics and Dynamics
J = robot.getBodyJacobian( q, nq );
J_l = J( 1:3, : );

q_sym = sym( 'q', [nq, 1 ] );
J_sym = robot.getBodyJacobian( q_sym, nq );
J_l_sym = J_sym( 1:3, : );
J_l_sym_T = J_l_sym';

% Partial derivative of J_T_sym wrt. q_sym
dJ_l_T_dq1 = diff( J_l_sym_T, q_sym(1) );
dJ_l_T_dq2 = diff( J_l_sym_T, q_sym(2) );
dJ_l_T_dq3 = diff( J_l_sym_T, q_sym(3) );
dJ_l_T_dq4 = diff( J_l_sym_T, q_sym(4) );
dJ_l_T_dq5 = diff( J_l_sym_T, q_sym(5) );
dJ_l_T_dq6 = diff( J_l_sym_T, q_sym(6) );
dJ_l_T_dq7 = diff( J_l_sym_T, q_sym(7) );

% Evaluate and store in matrix
dJ_l_T_dq = zeros( nq, 3, nq );
dJ_l_T_dq( :, :, 1 ) = double( subs( dJ_l_T_dq1, q_sym, q ) );
dJ_l_T_dq( :, :, 2 ) = double( subs( dJ_l_T_dq2, q_sym, q ) );
dJ_l_T_dq( :, :, 3 ) = double( subs( dJ_l_T_dq3, q_sym, q ) );
dJ_l_T_dq( :, :, 4 ) = double( subs( dJ_l_T_dq4, q_sym, q ) );
dJ_l_T_dq( :, :, 5 ) = double( subs( dJ_l_T_dq5, q_sym, q ) );
dJ_l_T_dq( :, :, 6 ) = double( subs( dJ_l_T_dq6, q_sym, q ) );
dJ_l_T_dq( :, :, 7 ) = double( subs( dJ_l_T_dq7, q_sym, q ) );

%% Stiffness matrices
K_l = zeros( 3, 3 );
K_l(1,1) = 402;
K_l(2,2) = 311;
K_l(3,3) = 167;

% K_r = zeros( 3, 3 );
% K_r(1,1) = 54;
% K_r(2,2) = 61;
% K_r(3,3) = 11;

%% External force and correction term K_Cs = \Gamma^k_{ij} * F_k
f = [ 15, 22, 14 ]';
% m = [ 5, 7, 11 ]';

% For moment part: K_m_CS = 0.5 * vec_to_so3( m );
K_f_CS = [ 0, -f(3), f(2); ...
           f(3), 0, -f(1); ...
            -f(2), f(1), 0 ];

% For 6D
% F = [ f; m ];
% K_CS = [ 0, 0, 0, 0, -F(3), F(2); ...
%         0, 0, 0, F(3), 0, -F(1); ...
%         0, 0, 0, -F(2), F(1), 0; ...
%         0, 0, 0, 0, -F(6)/2, F(5)/2; ...
%         0, 0, 0, F(6)/2, 0, -F(4)/2; ...
%         0, 0, 0, -F(5)/2, F(4)/2, 0 ];

%% Kinematic Stiffnesses (Chen, 2000, eq. 6) 
K_kin = [ dJ_l_T_dq( :, :, 1 ) * f, dJ_l_T_dq( :, :, 2 ) * f, dJ_l_T_dq( :, :, 3 ) * f, dJ_l_T_dq( :, :, 4 ) * f, ...
    dJ_l_T_dq( :, :, 5 ) * f, dJ_l_T_dq( :, :, 6 ) * f, dJ_l_T_dq( :, :, 7 ) * f ];

%% JS Stiffness Tensors
% Translations
Kq_l = J_l' * K_l * J_l;     

% Rotations
% Kq_r = J( 4:6, : )' * K_r * J( 4:6, : );

% Kinematic Stiffness to joint space
Kq_CS = J_l' * K_f_CS * J_l;

% Final joint space stiffness
Kq = Kq_l + K_kin                

%% SVD
[ U, Kq_bar, V ] = svd( Kq );
disp( Kq_bar )                      % Not singular






