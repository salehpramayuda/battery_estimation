function last_ekf(block)
% Level-2 MATLAB file S-Function 
% Discrete function for parameter estimation using
% extended kalman filter 
% estimate capacitor charge
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 6; %Delta, x0, P0, Z_k, W_k, Cg
  %% Register number of input and output ports
  block.NumInputPorts  = 3;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  %i_b(k)
  block.InputPort(1).Dimensions        = 1;
  block.InputPort(1).DirectFeedthrough = true;
  block.InputPort(1).SamplingMode='Sample'; 
  %i_l_k(k)
  block.InputPort(2).Dimensions        = 1;
  block.InputPort(2).DirectFeedthrough = true;
  block.InputPort(2).SamplingMode='Sample';   
  
  %y_k(k)
  block.InputPort(3).Dimensions        = 1;
  block.InputPort(3).DirectFeedthrough = true;
  block.InputPort(3).SamplingMode='Sample'; 
  
  
  %x(k)
  block.OutputPort(1).Dimensions       = [2,1];
  block.OutputPort(1).SamplingMode='Sample'; 
  
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
  block.NumDworks = 2;
  
  %states x
  block.Dwork(1).Name               = 'x'; 
  block.Dwork(1).Dimensions         = 1;
  block.Dwork(1).DatatypeID         = 0;
  block.Dwork(1).Complexity         = 'Real';
  block.Dwork(1).UsedAsDiscState    = true;
  
  %covariance matrix P
  block.Dwork(2).Name               = 'P';
  block.Dwork(2).Dimensions         = 1;
  block.Dwork(2).DatatypeID         = 0;
  block.Dwork(2).Complexity         = 'Real';
  block.Dwork(2).UsedAsDiscState    = true;

%endfunction

function InitConditions(block)

  %% Initialize Dwork
  %set x_prediction_k(0) and P_prediction_k(0) while calculating
  %x_k and P_k for the first iteration 
  x_p_k_0       = block.DialogPrm(2).data;
  P_p_k_0       = block.DialogPrm(3).data;

  block.Dwork(1).Data   = x_p_k_0;
  block.Dwork(2).Data   = P_p_k_0;
  
%endfunction

function Output(block)
  %assemble necessary elements
  Cg        = block.DialogPrm(6).data;
  Z_k       = block.DialogPrm(3).data;
  i_b_k     = block.InputPort(1).data;
  i_l_k     = block.InputPort(2).data;
  y_k       = block.InputPort(3).data;
  
  x_p_k     = block.Dwork(1).Data;
  P_p_k     = block.Dwork(2).Data;
  
% calculate kalman-gain and real covariance matrix

  y_hat = x_p_k/Cg;
  C = 1/Cg;
  
  K = P_p_k*C'/(Z_k + C*P_p_k*C');
  xk = x_p_k + K*(y_k - y_hat);
  Pk = P_p_k - K*C*P_p_k;
  
  %assign result to Dwork
  block.Dwork(1).Data = xk;
  block.Dwork(2).Data = Pk;
    
  %output estimated states  
  block.OutputPort(1).Data = [block.Dwork(1).Data; y_hat];
  
%endfunction

function Update(block)
  %update the states
  xk        = block.Dwork(1).Data;
  Pk        = block.Dwork(2).Data;
  Delta     = block.DialogPrm(1).Data;
  W_k       = block.DialogPrm(5).Data;
  i_b_k     = block.InputPort(1).Data;
  i_l_k     = block.InputPort(2).Data;
  
  %calculate prediction
  Ak = [1];
  
  x_p_kp1 = [xk + Delta*(i_l_k-i_b_k)];
  
  P_p_kp1 = Ak*Pk*Ak' + W_k;
 
  %update Dwork
  block.Dwork(1).Data   = x_p_kp1;
  block.Dwork(2).Data   = P_p_kp1;
 
%endfunction