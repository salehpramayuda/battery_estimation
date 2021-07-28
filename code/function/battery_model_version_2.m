function battery_model_version_2(block)
% Level-2 MATLAB file S-Function 
% Model Lead-Acid Battery based of ... et al (201x)
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 3; %x_conv, soc-uoc_polynome, Qe
  %% Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  %u(k)
  block.InputPort(1).Dimensions        = 1;
  block.InputPort(1).DirectFeedthrough = false;
  %block.InputPort(1).SamplingMode='Sample'; 

  %y(k)
  block.OutputPort(1).Dimensions       = 4;
  %block.OutputPort(1).SamplingMode='Sample'; 

  
  %% Set block sample time to inherited
  block.SampleTimes = [0 0];
  
  %% Setup Dwork
  block.NumContStates = 2;
  
  %% Set the block simStateCompliance to default 
  %% (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Derivatives',             @Derivatives);  
  
%endfunction


function InitConditions(block)

  %% Initialize Dwork
  block.ContStates.Data = block.DialogPrm(1).Data(4:5);
  
%endfunction

function Output(block)
  %assemble essential elements  
  param     = block.DialogPrm(1).Data(1:3);
  a         = block.DialogPrm(2).Data;
  x         = block.ContStates.Data;
  u         = block.InputPort(1).Data;
  
  u_oc = a(1)*x(1)^2+a(2)*x(1)+a(3);
  
  %calculate y(t)
  block.OutputPort(1).Data = [param(3)*(u-u_oc-x(2)), u_oc, x(1), x(2)];
  
%endfunction

function Derivatives(block)
  %calculate x_dot 
  u     = block.InputPort(1).Data;
  param = block.DialogPrm(1).Data(1:3);
  a     = block.DialogPrm(2).Data;
  Q     = block.DialogPrm(3).Data;
  x     = block.ContStates.Data;
  u_oc  = a(1)*x(1)^2+a(2)*x(1)+a(3);
  
  x_dashdot = [param(3)*(u-u_oc-x(2))/Q;
      -param(1)*(x(2)*param(2)-param(3)*(u-u_oc-x(2)))];
  block.Derivatives.Data = x_dashdot;
  
 %endfunction
