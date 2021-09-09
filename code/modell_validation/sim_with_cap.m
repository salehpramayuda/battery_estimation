function sim_with_cap(block)
% Level-2 MATLAB file S-Function 
% Discrete function for parameter estimation using extended kalman filter 
% Input: battery current; Output: voltage over battery (and parallel
% capacitor)
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 5; %Delta, x_conv, a, Qe, Cg
  %% Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 2;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  % u(k) battery current
  block.InputPort(1).Dimensions         = 1;
  block.InputPort(1).DirectFeedthrough  = true;
  block.InputPort(1).SamplingMode       ='Sample'; 
  
  % y(k) estimated voltage
  block.OutputPort(1).Dimensions        = 3;
  block.OutputPort(1).SamplingMode      ='Sample';
  
  % estimated currents
  block.OutputPort(2).Dimensions        = 2;
  block.OutputPort(2).SamplingMode      ='Sample';

  
  %% Set block sample time to inherited
  Delta = block.DialogPrm(1).Data;
  block.SampleTimes = [Delta 0];
  
  %% Set the block simStateCompliance to default 
  %% (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Update',                  @Update);  
%endfunction

function DoPostPropSetup(block)

  %% Setup Dwork
  block.NumDworks = 1;
  
  % States x
  block.Dwork(1).Name               = 'x'; 
  block.Dwork(1).Dimensions         = 3;
  block.Dwork(1).DatatypeID         = 0;
  block.Dwork(1).Complexity         = 'Real';
  block.Dwork(1).UsedAsDiscState    = true;
  
%endfunction

function InitConditions(block)

  %% Initialize Dwork

  block.Dwork(1).data   = block.DialogPrm(2).data(4:6);
%endfunction

function Output(block)
    %% Output the simulated voltage
    param   = block.DialogPrm(2).data(1:3);
    a       = block.DialogPrm(3).data;
    x       = block.Dwork(1).data;
    i_L     = block.InputPort(1).data;
    i_b     = x(3);
    
    block.OutputPort(1).data = x([1, 1, 2]);
    block.OutputPort(1).data(1) = i_b/param(3)+x(2)+a(1)*x(1)^2+a(2)*x(1)+a(3);
    block.OutputPort(2).data = [i_b; i_L-i_b];
%endfunction

function Update(block)
    %% Access important data
    i_L     = block.InputPort(1).Data;
    Delta   = block.DialogPrm(1).Data;
    param   = block.DialogPrm(2).Data(1:3);
    a       = block.DialogPrm(3).Data;
    Qe      = block.DialogPrm(4).Data;
    Cg      = block.DialogPrm(5).Data;
    
    x       = block.Dwork(1).Data;      % x = [SoC, u_rc, i_b, u_L] 
    %% Calculate next step for every state
    dsoc    = 1/Qe*x(3);
    durc    = param(1)*(x(3)-x(2)*param(2));
    dib     = param(3)*(i_L/Cg + x(2)*param(1)*param(2)-x(3)*(1/Cg+param(3)+...
            (2*a(1)*x(1)+a(2))/Qe));
    %dul     = 1/Cg*(i_L-x(3));
    
    block.Dwork(1).Data = [x(1)+Delta*dsoc; x(2)+Delta*durc;
            x(3)+Delta*dib];
    %block.Dwork(1).Data = [x(1)+Delta*dsoc; x(2)+Delta*durc;
            %x(3)+Delta*dib; x(4)+Delta*dul];


%endfunction
