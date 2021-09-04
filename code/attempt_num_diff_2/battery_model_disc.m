function battery_model_disc(block)
% Level-2 MATLAB file S-Function 
% Discrete function for parameter estimation using extended kalman filter 
% Input: battery current; Output: voltage over battery (and parallel
% capacitor)
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 4; %Delta, x_conv, y_km1, Qe, Cg
  %% Register number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  % u(k) battery current
  block.InputPort(1).Dimensions         = 1;
  block.InputPort(1).DirectFeedthrough  = true;
  block.InputPort(1).SamplingMode       ='Sample'; 
  
  % y_km1(k) measured voltage
  block.InputPort(2).Dimensions         = 1;
  block.InputPort(2).DirectFeedthrough  = true;
  block.InputPort(2).SamplingMode       ='Sample'; 
  
  % y(k) estimated voltage
  block.OutputPort(1).Dimensions        = 1;
  block.OutputPort(1).SamplingMode      ='Sample';

  
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
  % Set x_prediction_k(0) and P_prediction_k(0) while calculating
  % x_k and P_k for the first iteration 


  block.Dwork(1).data   = block.DialogPrm(2).data(5:7);
%endfunction

function Output(block)  

  % Output estimated states  
  block.OutputPort(1).data = block.Dwork(1).data(3);
  
%endfunction

function Update(block)
  % Access data for state update
  Delta = block.DialogPrm(1).data;
  param = block.DialogPrm(2).data(1:4);
  Qe    = block.DialogPrm(3).data;
  Cg    = block.DialogPrm(4).data;
  uk    = block.InputPort(1).data;
  y_km1 = block.InputPort(2).data;
  xk    = block.Dwork(1).data;
  i_c   = Cg*(xk(3)-y_km1)/Delta;
  
  % Calculate next state
  soc = xk(1) + Delta*uk/Qe;
  urc = xk(2) + Delta*param(1)*(uk-xk(2)*param(2));
  u_L = xk(3) + Delta*i_c/Cg;
 
  % Update Dwork
  block.Dwork(1).data   = [soc; urc; u_L];
%endfunction
