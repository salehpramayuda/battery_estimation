function sim_wo_cap(block)
% Level-2 MATLAB file S-Function 
% Discrete function for parameter estimation using extended kalman filter 
% Input: battery current; Output: voltage over battery (and parallel
% capacitor)
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 3; %x_conv, a, Qe
  %% Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  % u(k) battery current
  block.InputPort(1).Dimensions         = 1;
  block.InputPort(1).DirectFeedthrough  = false;
  %block.InputPort(1).SamplingMode       ='Sample'; 
  
  % y(k) estimated voltage
  block.OutputPort(1).Dimensions        = 3;
  %block.OutputPort(1).SamplingMode      ='Sample';

  
  %% Set block sample time to inherited
  block.SampleTimes = [0 0];
  block.NumContStates = 2;
  
  %% Set the block simStateCompliance to default 
  %% (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Derivatives',             @Derivative);  
%endfunction

function InitConditions(block)

  %% Initialize Dwork

  block.ContStates.data   = block.DialogPrm(1).data(4:5);
%endfunction

function Output(block)
    %% Output the simulated voltage
    u       = block.InputPort(1).data;
    param   = block.DialogPrm(1).data(1:3);
    a       = block.DialogPrm(2).data;
    x       = block.ContStates.data;
    
    u_L     = u/param(3) + x(2) + a(1)*x(1)^2+a(2)*x(1)+a(3);
    
    block.OutputPort(1).Data = [u_L; x(1); x(2)];
%endfunction

function Derivative(block)
    %% Access important data
    u       = block.InputPort(1).data;
    param   = block.DialogPrm(1).data(1:3);
    Qe      = block.DialogPrm(3).data;
    
    x       = block.ContStates.data;      % x = [SoC, u_rc]
    
    %% Calculate next step for every state
    dsoc    = u/Qe;
    durc    = param(1)*(u-x(2)*param(2));
    
    block.Derivatives.data = [dsoc; durc];
%endfunction
