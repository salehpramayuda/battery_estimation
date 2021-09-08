function simulate_with_cap(block)
    % time continuous level 2 s-function to simulate
    % "open-circuit" voltage of battery with an external
    % capacitor (see accurate model)
    % 3 state (u_L, u_rc, soc)
    setup(block);
    
%endfunction

function setup(block)
    
    block.NumDialogPrms = 4; %x0 (initial state + param), a, Qe, Cg
    %% Register number of input and output ports
    block.NumInputPorts  = 1;
    block.NumOutputPorts = 2;

    %% Setup functional port properties to dynamically
    %  inherited.
    block.SetPreCompInpPortInfoToDynamic;
    block.SetPreCompOutPortInfoToDynamic;
  
    % u battery current
    block.InputPort(1).Dimensions        = 1;
    block.InputPort(1).DirectFeedthrough = true;
    block.InputPort(1).SamplingMode     ='Sample'; 
  
    % y_hat simulated current
    block.OutputPort(1).Dimensions      = 1;
    block.OutputPort(1).SamplingMode    ='Sample';
    
    % currents
    block.OutputPort(2).Dimensions      = [2,1];
    block.OutputPort(2).SamplingMode    ='Sample'; 
  
    %% Set block sample time to continuous
    block.SampleTimes = [0 0];
    
    %% Setup Dwork
    block.NumContStates = 4;
  
    %% Set the block simStateCompliance to default 
    %% (i.e., same as a built-in block)
    block.SimStateCompliance = 'DefaultSimState';

    %% Register methods
    block.RegBlockMethod('InitializeConditions',    @InitConditions);  
    block.RegBlockMethod('Outputs',                 @Output);  
    block.RegBlockMethod('Derivatives',              @Derivative);  
    
%endfunction

function InitConditions(block)
    %% Initialize continuous states SoC, u_rc, i_b, u_L
    %block.ContStates.Data = [block.DialogPrm(1).Data(6); block.DialogPrm(1).Data(4:5)];
    block.ContStates.Data = block.DialogPrm(1).Data(4:7);
%endfunction

function Output(block)
    %% Output the simulated voltage
    i_L     = block.InputPort(1).Data;
    i_b     = block.ContStates.Data(3);
    
    block.OutputPort(1).Data = block.ContStates.Data(4);
    block.OutputPort(2).Data = [i_b, i_L-i_c];
%endfunction

function Derivative(block)
    %% Access important data
    i_L     = block.InputPort(1).Data;
    param   = block.DialogPrm(1).Data(1:3);
    a       = block.DialogPrm(2).Data;
    Qe      = block.DialogPrm(3).Data;
    Cg      = block.DialogPrm(4).Data;
    
    x       = block.ContStates.Data;      % x = [u_L, SoC, u_rc] 
    %% Calculate derivatives for every state
    dsoc    = 1/Qe*x(3);
    durc    = param(1)*(x(3)-x(2)*param(2));
    dib     = param(3)*(i_L/Cg + x(2)*param(1)*param(2)-x(3)*(1/Cg+param(3)+...
            (2*a(1)*x(1)+a(2))/Qe));
    dul     = 1/Cg*(i_L-x(3));
    
    block.Derivatives.Data = [dsoc;durc;dib;dul];

%endfunction
