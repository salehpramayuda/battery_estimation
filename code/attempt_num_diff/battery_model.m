function battery_model(block)
% Level-2 MATLAB file S-Function 
% Model Lead-Acid Battery based of ... et al (201x)
% modified to be suited for Vrontos (2020)
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 4; %x_conv, soc-uoc_polynome, Qe, Cg
  %% Register number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  % u(k)
  block.InputPort(1).Dimensions         = 1;
  block.InputPort(1).DirectFeedthrough  = true;
  %block.InputPort(1).SamplingMode='Sample';
  
  % y(k-1)
  block.InputPort(2).Dimensions         = 1;
  block.InputPort(2).DirectFeedthrough  = true;

  % y(k)
  block.OutputPort(1).Dimensions        = 4;
  %block.OutputPort(1).SamplingMode='Sample'; 

  
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

  %% Initialize Dwork
  block.ContStates.Data = block.DialogPrm(1).data(5:7);
  
%endfunction

function Output(block)
  % Assemble essential elements  
  x         = block.ContStates.data;
  
  % Output y
  block.OutputPort(1).data = x(3);
  
%endfunction

function Derivatives(block)
  %calculate x_dot 
  u     = block.InputPort(1).data;
  y_km1 = block.InputPort(2).data;
  param = block.DialogPrm(1).data(1:4);
  a     = block.DialogPrm(2).data;
  Qe    = block.DialogPrm(3).data;
  Cg    = block.DialogPrm(4).data;
  x     = block.ContStates.data;
  
  i_c   = Cg*(x(3)-y_km1)/
  
  x_dashdot = [param(3)*(u-u_oc-x(2))/Q; -param(1)*(x(2)*param(2)-param(3)*(u-u_oc-x(2)))];
  block.Derivatives.data = x_dashdot;
  
 %endfunction
