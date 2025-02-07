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
% anim = Animation( 'Dimension', 3, 'xLim', [-0.7,0.7], 'yLim', [-0.7,0.7], 'zLim', [0,1.4] );
% anim.init( );
% anim.attachRobot( robot )  

%% Update kinematics
robot.updateKinematics( robot.q_init );
% anim.update( 0 );

%% Kinematics and Dynamics
q = robot.q_init;
J = robot.getHybridJacobian( q );
J = J( 4:6, : );

q_sym = sym( 'q', [nq, 1 ] );
J_sym = robot.getHybridJacobian( q_sym );
J_sym = J_sym( 4:6, : );
J_T_sym = J_sym';

% Partial derivative of J_T_sym wrt. q_sym
dJ_T_dq1 = diff( J_T_sym, q_sym(1) );
dJ_T_dq2 = diff( J_T_sym, q_sym(2) );
dJ_T_dq3 = diff( J_T_sym, q_sym(3) );
dJ_T_dq4 = diff( J_T_sym, q_sym(4) );
dJ_T_dq5 = diff( J_T_sym, q_sym(5) );
dJ_T_dq6 = diff( J_T_sym, q_sym(6) );
dJ_T_dq7 = diff( J_T_sym, q_sym(7) );

% Evaluate and store in matrix
dJ_T_dq = zeros( nq, 3, nq );
dJ_T_dq( :, :, 1 ) = double( subs( dJ_T_dq1, q_sym, q ) );
dJ_T_dq( :, :, 2 ) = double( subs( dJ_T_dq2, q_sym, q ) );
dJ_T_dq( :, :, 3 ) = double( subs( dJ_T_dq3, q_sym, q ) );
dJ_T_dq( :, :, 4 ) = double( subs( dJ_T_dq4, q_sym, q ) );
dJ_T_dq( :, :, 5 ) = double( subs( dJ_T_dq5, q_sym, q ) );
dJ_T_dq( :, :, 6 ) = double( subs( dJ_T_dq6, q_sym, q ) );
dJ_T_dq( :, :, 7 ) = double( subs( dJ_T_dq7, q_sym, q ) );

%% Stiffness matrix
K = zeros( 3, 3 );
K(1,1) = 40;
K(2,2) = 31;
K(3,3) = 16;

%% External force and correction term C = \Gamma^k_{ij} * F_k
m = [ 1.5, 2.2, 1.4 ]';

% No correction terms!

% For 6D
% C = [ 0, 0, 0, 0, -F(3), F(2); ...
%         0, 0, 0, F(3), 0, -F(1); ...
%         0, 0, 0, -F(2), F(1), 0; ...
%         0, 0, 0, 0, -F(6)/2, F(5)/2; ...
%         0, 0, 0, F(6)/2, 0, -F(4)/2; ...
%         0, 0, 0, -F(5)/2, F(4)/2, 0 ];

% Chen 2004!!!!!!!!!!!!!!!!!!!!
K_CS = [ 0, m(3)/2, -m(2)/2; ...
       -m(3)/2, 0, m(1)/2; ...
       m(2)/2, -m(1)/2, 0 ];

%% Kinematic Stiffness (Chen, 2000, eq. 6) -> Symmetric for 3D!
K_kin = [ dJ_T_dq( :, :, 1 ) * m, dJ_T_dq( :, :, 2 ) * m, dJ_T_dq( :, :, 3 ) * m, dJ_T_dq( :, :, 4 ) * m, ...
    dJ_T_dq( :, :, 5 ) * m, dJ_T_dq( :, :, 6 ) * m, dJ_T_dq( :, :, 7 ) * m ];

%% JS Stiffness Tensor
K_js = J' * K * J + J' * K_CS * J + K_kin

[U, K_bar, V] = svd( K_js );




