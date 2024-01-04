% TITLE: Automated Airfoil Optimization Using XFoil V1.2.2
% AUTHOR: Andrew Smith
% UPDATED: 11/14/2023 
%
% INSTRUCTIONS: 
% 1. Download the XFoil Matlab interface from the site below:
%    https://www.mathworks.com/matlabcentral/fileexchange/49706-xfoil-interface-updated
% 2. Save this file into the same folder as the XFoil program and related 
%    files.  
% 3. Enter all necessary input parameters for the airfoil such as file save 
%    name, fluid parameters, angles of attack, etc (lines 34-57) .  
% 4. Call the "AUTO_XFOIL_xfoil" function using the initial airfoil input 
%    shape parameters shown below:
%    - A through D = shape parameters to tune (FLOAT)
%    - CAM = maximum camber percentage (%)
%    - LOC = location of maximum camber (%)
%    Example: AUTO_XFOIL(10,2,2,10,0.05,0.4)
% 5. Allow the program time to complete all simulations needed. Stop the
%    program at any point using CTRL+C. 
% 6. Once an optimal airfoil is found, repeat the analysis by decreasing
%    the step sizes for each parameter to tune the individual values if
%    desired. 
% 7. Perform an XFoil simulation on the final airfoil to validate the 
%    results found. 
%
% NOTES:
% - Add code to allow user to determine operating parameter to maximize
%   from XFOIL (e.g., CL, CL/CD, -CD). 

% ===================================================================================================================


function cl_avg = AUTO_XFOIL(A, B, C, D, CAM, LOC)
% INPUT PARAMETERS - XFOIL SIMULATION
% File Name & Location
fileSaveName = 'C:\Users\andre\XFOIL6.99\ARS_FOIL.dat'; % location where .DAT airfoil file will save. *** note: file will be overwritten with each simulation. ***
% Angle of Attack
aoa_start = 2; 
aoa_end   = 7; % for single angle of attack, set aoa_end = aoa_start
aoa_step  = 1;
% Fluid Parameters
reynolds = 1e6; 
mach     = 0  ; % set 0 for incompressible
% Simulation Parameters
target = "CL"; % ["CL", "CL/CD", or "-CD"] parameter to maximize from XFoil
iter = 500; % maximum XFoil iterations

% INPUT PARAMETERS - AUTOMATION VARIABLES
% Maximum Iterations
iter_max = 10;
% Parameters Step Size
A_step = 1;
B_step = 1;
C_step = 1;
D_step = 1;
CAM_step = 0.001;
LOC_step = 0.01;

% Display input values
disp("Shape Parameters: A = " + A+ ", B = " + B + ", C = " + C + ", D = " + D + ", CAM = " + CAM + ", LOC = " + LOC);
disp("Angle of Attack: alpha = " + aoa_start + " : " + aoa_step + " : " + aoa_end);
disp(" ");
disp("Generating airfoil...");
% create airfoil based on parameters
my_shape.A = A;
my_shape.B = B;
my_shape.C = C;
my_shape.D = D;
my_camber.max = CAM;
my_camber.loc = LOC;
[x,y] = makeFoil(my_shape,my_camber,fileSaveName); 

disp("Airfoil successfully generated.");
cont_1 = input("Continue? [1 = YES, 0 = NO] : ");
if cont_1 == 0
    disp("Program terminated.");
    return % Exit function
end

% Find Initial target (e.g. CL, CL/CD).  
disp(" ");
disp("Determining Current target value...");

CLs_current = []; % CL used to refer to target value
for AOA = aoa_start:aoa_step:aoa_end % Define range of angles of attack
    if target == "CL"
        CLs_current = [CLs_current, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
    elseif target == "CL/CD"
        CLs_current = [CLs_current, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
    elseif target == "-CD"
        CLs_current = [CLs_current, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
    end
end
cl_initial = mean(CLs_current, "all"); % Average target value
if target == "CL"
    disp("Current CL : " + cl_initial);
elseif target == "CL/CD"
    disp("Current CL/CD : " + cl_initial);
elseif target == "-CD"
    disp("Current CD : " + (-1*cl_initial) );
end

% HILL CLIMBING ALGORITHM 
% determine current best parameters and output values 
CL_max = cl_initial; % current maximum (use initial airfoil to start)
A_max = A; % shape parameters at maximum value above (set at input parameters initially)
B_max = B;
C_max = C;
D_max = D;
CAM_max = CAM;
LOC_max = LOC;

num_simulations = iter_max * 12 * (ceil((aoa_end - aoa_start)/aoa_step) + 1); % maximum number of XFoil simulation needed based on inputs
disp(" ");
disp("Expected Number of XFoil Simulations: " + num_simulations );
cont_2 = input("Continue? [1 = YES, 0 = NO] : ");
if cont_2 == 0
    disp("Program terminated.");
    return % Exit function
end
disp("Optimization running. Press CTRL+C to stop at any point.")

i = 0; % current iteration
while i < iter_max  
    cl_current = CL_max; % set the current target value as the previous max (used to stop program if nothing increases CL_max from cl_current)
    
    % ================ Alter A parameter (-) ================
    A_new = A_max - A_step; 
    
    my_shape.A = A_new; % change A
    my_shape.B = B_max;
    my_shape.C = C_max;
    my_shape.D = D_max;
    my_camber.max = CAM_max;
    my_camber.loc = LOC_max;
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        A_max = A_new; % if better, change A
        CL_max = cl_new;
    end
        
    
    % ================ Alter A parameter (+) ================
    A_new = A_max + A_step;
    if CL_max == cl_new % If previous design was better, above code will retest current design parameters
        A_new = A_max + 2*A_step;
    end
    
    my_shape.A = A_new; % change A
    my_shape.B = B_max;
    my_shape.C = C_max;
    my_shape.D = D_max;
    my_camber.max = CAM_max;
    my_camber.loc = LOC_max;
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        A_max = A_new; % if better, change A
        CL_max = cl_new;
    end
    
    
    % ================ Alter B parameter (-) ================
    B_new = B_max - B_step; 
    
    my_shape.A = A_max; 
    my_shape.B = B_new; % change B
    my_shape.C = C_max;
    my_shape.D = D_max;
    my_camber.max = CAM_max;
    my_camber.loc = LOC_max;
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        B_max = B_new; % if better, change B
        CL_max = cl_new;
    end
        
    
    % ================ Alter B parameter (+) ================
    B_new = B_max + B_step; 
    if CL_max == cl_new % If previous design was better, above code will retest current design parameters
        B_new = B_max + 2*B_step;
    end
    
    my_shape.A = A_max; 
    my_shape.B = B_new; % change B
    my_shape.C = C_max;
    my_shape.D = D_max;
    my_camber.max = CAM_max;
    my_camber.loc = LOC_max;
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        B_max = B_new; % if better, change B
        CL_max = cl_new;
    end
    
    
    % ================ Alter C parameter (-) ================
    C_new = C_max - C_step; 
    
    my_shape.A = A_max; 
    my_shape.B = B_max; 
    my_shape.C = C_new; % change C
    my_shape.D = D_max;
    my_camber.max = CAM_max;
    my_camber.loc = LOC_max;
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        C_max = C_new; % if better, change C
        CL_max = cl_new;
    end
        
    
    % ================ Alter C parameter (+) ================
    C_new = C_max + C_step;
    if CL_max == cl_new % If previous design was better, above code will retest current design parameters
        C_new = C_max + 2*C_step;
    end
    
    my_shape.A = A_max; 
    my_shape.B = B_max; 
    my_shape.C = C_new; % change C
    my_shape.D = D_max;
    my_camber.max = CAM_max;
    my_camber.loc = LOC_max;
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        C_max = C_new; % if better, change C
        CL_max = cl_new;
    end
    
    
    % ================ Alter D parameter (-) ================
    D_new = D_max - D_step; 
    
    my_shape.A = A_max; 
    my_shape.B = B_max; 
    my_shape.C = C_max; 
    my_shape.D = D_new; % change D
    my_camber.max = CAM_max;
    my_camber.loc = LOC_max;
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        D_max = D_new; % if better, change D
        CL_max = cl_new;
    end
    
    
    % ================ Alter D parameter (+) ================
    D_new = D_max + D_step;
    if CL_max == cl_new % If previous design was better, above code will retest current design parameters
        D_new = D_max + 2*D_step;
    end
    
    my_shape.A = A_max; 
    my_shape.B = B_max; 
    my_shape.C = C_max; 
    my_shape.D = D_new; % change D
    my_camber.max = CAM_max;
    my_camber.loc = LOC_max;
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        D_max = D_new; % if better, change D
        CL_max = cl_new;
    end
    
    
    % ================ Alter CAM parameter (-) ================
    CAM_new = CAM_max - CAM_step; 
    
    my_shape.A = A_max; 
    my_shape.B = B_max; 
    my_shape.C = C_max; 
    my_shape.D = D_max; 
    my_camber.max = CAM_new; % change CAM
    my_camber.loc = LOC_max;
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        CAM_max = CAM_new; % if better, change CAM
        CL_max = cl_new;
    end
    
    
    % ================ Alter CAM parameter (+) ================
    CAM_new = CAM_max + CAM_step;
    if CL_max == cl_new % If previous design was better, above code will retest current design parameters
        CAM_new = CAM_max + 2*CAM_step;
    end
    
    my_shape.A = A_max; 
    my_shape.B = B_max; 
    my_shape.C = C_max; 
    my_shape.D = D_max; 
    my_camber.max = CAM_new; % change CAM
    my_camber.loc = LOC_max;
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        CAM_max = CAM_new; % if better, change CAM
        CL_max = cl_new;
    end
    
    
    % ================ Alter LOC parameter (-) ================
    LOC_new = LOC_max - LOC_step; 
    
    my_shape.A = A_max; 
    my_shape.B = B_max; 
    my_shape.C = C_max; 
    my_shape.D = D_max; 
    my_camber.max = CAM_max; 
    my_camber.loc = LOC_new; % change LOC
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        LOC_max = LOC_new; % if better, change LOC
        CL_max = cl_new;
    end
    
    
    % ================ Alter LOC parameter (+) ================
    LOC_new = LOC_max + LOC_step;
    if CL_max == cl_new % If previous design was better, above code will retest current design parameters
        LOC_new = LOC_max + 2*LOC_step;
    end
    
    my_shape.A = A_max; 
    my_shape.B = B_max; 
    my_shape.C = C_max; 
    my_shape.D = D_max; 
    my_camber.max = CAM_max; 
    my_camber.loc = LOC_new; % change LOC
    [x,y] = makeFoil(my_shape,my_camber,fileSaveName);
    
    CLs_new = [];
    for AOA = aoa_start:aoa_step:aoa_end
        if target == "CL"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL) ]; % Add CL to list of values
        elseif target == "CL/CD"
           CLs_new = [CLs_new, (xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CL)/(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add CL/CD to list of values
        elseif target == "-CD"
           CLs_new = [CLs_new, -(xfoil('ARS_FOIL.dat', AOA, reynolds, 0, 'pane', 'oper iter ' + string(iter) ).CD) ]; % Add -CD to list of values
        end
    end
    cl_new = mean(CLs_new, "all");
    
    if cl_new > CL_max
        LOC_max = LOC_new; % if better, change LOC
        CL_max = cl_new;
    end
    
    
    % ================ Display Results ================
    if CL_max == cl_current % if CL_max has not changed from current, then current is local maximum
        break
    end
    
    disp(" ");
    disp("Iteration " + (i+1) + " of " + iter_max);
    disp("A = " + A_max + ", B = " + B_max + ", C = " + C_max + ", D = " + D_max + ", CAM = " + CAM_max + ", LOC = " + LOC_max);
    if target == "CL"
        disp("CL = " + CL_max);
    elseif target == "CL/CD"
        disp("CL/CD = " + CL_max);
    else
        disp("CD = " + (-1*CL_max) );
    end
    i = i + 1; % increase iterations (used to halt program)
end

disp(" ");
disp("Best parameters found: ");
disp("A = " + A_max);
disp("B = " + B_max);
disp("C = " + C_max);
disp("D = " + D_max);
disp("CAM = " + CAM_max);
disp("LOC = " + LOC_max);

disp(" ");
disp("Optimized target value: ");
disp("CL = " + CL_max);

disp(" ");
disp("To rerun the optimization, copy the function call below: ");
disp("AUTO_XFOIL("+A_max+","+B_max+","+C_max+","+D_max+","+CAM_max+","+LOC_max+")");



% ===================================================================================================================



function [x,y] = makeFoil(shape,camber,fileSaveName)

%% DEFINE FOIL SHAPE
load('foilLibrary.mat','A_mean','U_a','S_a','V_a','sampleX')

[~,xloc] = min(abs(sampleX-camber.loc)); %nearest integer where camber pieces meet
camber.y1 = (camber.max/camber.loc^2) * (2*camber.loc*sampleX(1:xloc) - sampleX(1:xloc).^2);                        %camber piece 1 (front)
camber.y2 = camber.max/(1-camber.loc)^2*((1-2*camber.loc)+2*camber.loc*sampleX(xloc+1:end)-sampleX(xloc+1:end).^2); %camber piece 2 (back)
camber.y = [camber.y1 camber.y2];

coefficients = [shape.A shape.B shape.C shape.D].*[0.0008 0.0023 -0.0126 -0.0038]; 

%% ENSURE FOIL IS PHYSICAL (SURFACE CAN'T CROSS CHORDLINE)
U_b = U_a*S_a;
U = U_b(:,1:length(coefficients));
S = S_a(1:length(coefficients),1:length(coefficients));
tmin = 0.000;
options = optimset('Display','off');
% coefficientsN = quadprog(S'*S,-coefficients*S'*S,-U*S,A_mean,[],[],[],[],[],options)';
options = optimset('Display','off');
coefficientsN = quadprog(U'*U,-U'*U*coefficients',-U,A_mean - tmin,[],[],[],[],[],options);

modeNum = length(coefficientsN);
    
%% MAKE FOIL
foil_SVD = zeros(length(sampleX),1);
for k = 1:modeNum
    mode = U_a(:,k)*S_a(k,k)*coefficientsN(k)';
    mode(1) = 0;
    mode(end) = 0;
    foil_SVD = foil_SVD + mode;
end
A_mean(1) = 0;
A_mean(end) = 0;
foil_SVD = foil_SVD + A_mean;

%% PLOT FOIL
plot(sampleX,foil_SVD'+camber.y,'.-r')
hold on
plot(sampleX,-foil_SVD'+camber.y,'.-r')
set(gca,'fontname','Source Sans Pro Light','fontsize',12,'ygrid','on','xgrid','on')
xlim([0 1])
set(gcf,'position',[520   600   560   198])
plot(sampleX,camber.y,'--k')
axis equal

%% SAVE FOIL
x = [sampleX(end:-1:1) sampleX];
top = foil_SVD'+camber.y;
bottom = -foil_SVD'+camber.y;
y = [top(end:-1:1) bottom];
mat = [x' y'];
save(fileSaveName,'mat','-ascii')