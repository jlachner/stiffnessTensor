% [Project]        Robot Simulator - Snake
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
simTime = 0.01;        % Total simulation time
t  = 0;             % The current time of simulation
dt = 0.01;          % Time-step of simulation

%% Initialize the robot

% Geometric and Inertial Parameters of SnakeBot
nq = 3;         % The number of linkages of the Snakebot
m  = 1;         % The   mass of the each link
l  = 1;         % The length of the each link

m_arr = m * ones( 1, nq );  % The mass   array to construct SnakeBot
l_arr = l * ones( 1, nq );  % The length array to construct SnakeBot

% Construct a 5-DOF SnakeBot
robot = SnakeBot( nq, m_arr, l_arr );
robot.init( )

% Attach the 5-DOF SnakeBot to animation for visualization
anim = Animation( 'Dimension', 2, 'xLim', [-1.5,6.5], 'yLim', [-4,4], 'isSaveVideo', false, 'VideoSpeed', 0.5 );
anim.init( )
anim.attachRobot( robot )


%% Initialization of Animation

% DO NOT CHANGE
% Changing the degrees to radian
q = deg2rad( [ 30, 22, -15 ]' );
q_sym = sym( 'q', [ nq, 1 ] );
dq = zeros( nq, 1 );
% Update robot kinematics with q_deg array
% Also get the end-effector's H matrix
robot.updateKinematics( q );
H_EE = robot.getForwardKinematics( q );
H_ini = H_EE;

% Update animation
anim.update( 0 );


%% Running the main-loop of simulation

while t <= simTime


    % Get the p array from the SE(3) Matrix
    p_EE = H_EE( 1:3, 4 );

    % Get Hybrid Jacobian
    J = robot.getHybridJacobian( q );
    J_sym = robot.getHybridJacobian( q_sym );
    J_T_sym = J_sym';

    % Get Mass Matrix
    M = robot.getMassMatrix( q );

    % Stiffness matrix
    K_mat = zeros( 6, 6 );
    K_mat( 1, 1 ) = 15;
    K_mat( 2, 2 ) = 11;
    K_mat( 3, 3 ) = 27;

    % Curvature correction term C = \Gamma^k_{ij} * F_k
    F = [ 1.5, 2.2, 0, 0, 0, 0 ]';

    % NOT NEEDED FOR 3D
    C = [ 0, 0, 0, 0, -F(3), F(2); ...
        0, 0, 0, F(3), 0, -F(1); ...
        0, 0, 0, -F(2), F(1), 0; ...
        0, 0, 0, 0, -F(6)/2, F(5)/2; ...
        0, 0, 0, F(6)/2, 0, -F(4)/2; ...
        0, 0, 0, -F(5)/2, F(4)/2, 0 ];

    % Cart. stiffness tensor
    K = K_mat - C;                      % NOT NEEDED FOR 3D SINCE LIE BRACKETS VANISH!
    K = K_mat;

    % Partial derivative of Jacobian
    dJ_T_dq1 = diff( J_T_sym, q_sym( 1 ) );
    dJ_T_dq2 = diff( J_T_sym, q_sym( 2 ) );
    dJ_T_dq3 = diff( J_T_sym, q_sym( 3 ) );

    % Store in one matrix
    dJT_T_dq_mat = zeros( 3, 6, nq );
    dJT_T_dq_mat( :, :, 1) = double( subs( dJ_T_dq1, q_sym, q ) );
    dJT_T_dq_mat( :, :, 2) = double( subs( dJ_T_dq2, q_sym, q ) );
    dJT_T_dq_mat( :, :, 3) = double( subs( dJ_T_dq3, q_sym, q ) );

    K_kin = [ dJT_T_dq_mat( :, :, 1 ) * F, dJT_T_dq_mat( :, :, 2 ) * F, dJT_T_dq_mat( :, :, 3 ) * F ] 

    % JS stiffness tensor
    K_js = J' * K * J + K_kin

    % Superposition of torques
    tau = zeros( nq, 1 );

    % proceed one simulation step
    % rhs = M\(torque -c - g);        % if coriolis and gravity should be calculated
    rhs = M\(tau);                    % forgo coriolis and gravity for simulation


    [ q1, dq1 ] = func_symplecticEuler( q, dq, rhs, dt);
    q  =  q1;
    dq = dq1;


    % Update the linkage plot
    robot.updateKinematics( q );

    anim.update( t );


    % Get the forward kinematics of the EE
    H_EE = robot.getForwardKinematics( q );
    t = t + dt;

end
