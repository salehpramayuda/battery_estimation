function battery_model(block)
% Level-2 MATLAB S-Function
% Describes the battery model for validation purposes

setup(block)
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 3; %x_conv, Qe, Cg
  %% Register number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  %i_b(k) battery currrent
  block.InputPort(1).Dimensions        = 1;
  block.InputPort(1).DirectFeedthrough = true;
  block.InputPort(1).SamplingMode='Sample'; 
  
  %i_L(k) terminal current
  block.InputPort(2).Dimensions        = 1;
  block.InputPort(2).DirectFeedthrough = true;
  block.InputPort(2).SamplingMode='Sample'; 
  
  %x(k)
  block.OutputPort(1).Dimensions       = [3,1];
  block.OutputPort(1).SamplingMode='Sample'; 
  
  %% Set block sample time to inherited
  block.SampleTimes = [0 0];
  
   %% Setup Dwork
  block.NumContStates = 3;
  
  %% Set the block simStateCompliance to default 
  %% (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Derivatives',             @Derivatives);  
  
%endfunction


function InitConditions(block)
  block.ContStates.Data = block.DialogPrm(1).Data([4:5,7]);
  
%endfunction

function Output(block)
  %assemble essential elements
  Cg        = block.DialogPrm(2).Data;
  x         = block.ContStates.Data;
 
  %calculate y(t)
  block.OutputPort(1).Data = [x(3)/Cg; x(1); x(2)];
  
%endfunction

function Derivatives(block)
  %calculate x_dot 
  i_b   = block.InputPort(1).Data;
  i_l   = block.InputPort(2).Data;
  param = block.DialogPrm(1).Data(1:3);
  Q     = block.DialogPrm(3).Data;
  x     = block.ContStates.Data;
  
  x_dashdot = [param(1)*(i_b-x(1)*param(2)); i_b/Q; i_l-i_b];
  block.Derivatives.Data = x_dashdot;
  
 %endfunction
