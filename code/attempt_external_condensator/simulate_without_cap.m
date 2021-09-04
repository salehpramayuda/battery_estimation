function simulate_without_cap(block)
    % time continuous level 2 s-function to simulate
    % "open-circuit" voltage of battery without an external
    % capacitor (see accurate model)
    %  state (u_rc, soc)
    setup(block);
    
%endfunction

function setup(block)
    
    block.NumDialogPrms = 3; %x0 (initial state + param), a, Qe
    %% Register number of input and output ports
    block.NumInputPorts  = 1;
    block.NumOutputPorts = 1;

    %% Setup functional port properties to dynamically
    %  inherited.
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;
  
    % u battery current
    block.InputPort(1).Dimensions        = 1;
    block.InputPort(1).DirectFeedthrough = true;
    block.InputPort(1).SamplingMode='Sample'; 
  
    % y_hat simulated current
    block.OutputPort(1).Dimensions       = 1;
    block.OutputPort(1).SamplingMode='Sample'; 
  
    %% Set block sample time to continuous
    block.SampleTimes = [0 0];
    
    %% Setup Dwork
    block.NumContStates = 2;
  
    %% Set the block simStateCompliance to default 
    %% (i.e., same as a built-in block)
    block.SimStateCompliance = 'DefaultSimState';

    %% Register methods
    block.RegBlockMethod('InitializeConditions',    @InitConditions);  
    block.RegBlockMethod('Outputs',                 @Output);  
    block.RegBlockMethod('Derivatives',              @Derivative);  
    
%endfunction

function InitConditions(block)
    %% Initialize continuous states u_L, SoC, u_rc
    block.ContStates.Data = block.DialogPrm(1).Data(4:5);

%endfunction

function Output(block)
    %% Output the simulated voltage
    param   = block.DialogPrm(1).Data(1:3);
    a       = block.DialogPrm(2).Data;
    x       = block.ContStates.Data;
    i_b     = block.InputPort(1).Data;
    
    u_oc    = a(1)*x(1)^2+a(2)*x(1)+a(3);
    u_L     = u_oc + x(2) + i_b/param(3);
    
    block.OutputPort(1).Data = u_L;
%endfunction

function Derivative(block)
    %% Access important data
    i_b     = block.InputPort(1).Data;
    param   = block.DialogPrm(1).Data(1:3);
    a       = block.DialogPrm(2).Data;
    Qe      = block.DialogPrm(3).Data;
    x       = block.ContStates.Data;      % x = [SoC, u_rc]
    
    %% Calculate derivatives for every state
    dsoc    = 1/Qe*i_b;
    durc    = param(1)*(i_b-x(2)*param(2));
    
    block.Derivatives.Data = [dsoc;durc];

%endfunction
