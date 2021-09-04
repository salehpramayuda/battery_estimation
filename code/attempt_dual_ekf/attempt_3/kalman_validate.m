function kalman_validate(block)
% Level-2 MATLAB file S-Function 
% Discrete function for parameter estimation using
% extended kalman filter 
% input: battery current, output: voltage over battery and Capacitor of B-BC 
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 9; %Delta, x_0, P_0, soc-uoc_polynome, 
                        %Z_k, W_k, Qe, param
  %% Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  %u(k) battery currrent
  block.InputPort(1).Dimensions        = 1;
  block.InputPort(1).DirectFeedthrough = true;
  block.InputPort(1).SamplingMode='Sample'; 
  
%   %urc RC-voltage
%   block.InputPort(2).Dimensions        = 1;
%   block.InputPort(2).DirectFeedthrough = true;
%   block.InputPort(2).SamplingMode='Sample'; 
%   
%   %urc battery soc
%   block.InputPort(3).Dimensions        = 1;
%   block.InputPort(3).DirectFeedthrough = true;
%   block.InputPort(3).SamplingMode='Sample'; 
  
  %x(k)
  block.OutputPort(1).Dimensions       = 1;
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
  block.NumDworks = 3;
  
  %states x
  block.Dwork(1).Name               = 'x'; 
  block.Dwork(1).Dimensions         = 2;
  block.Dwork(1).DatatypeID         = 0;
  block.Dwork(1).Complexity         = 'Real';
  block.Dwork(1).UsedAsDiscState    = true;
  
  %covariance matrix P
  block.Dwork(2).Name               = 'P';
  block.Dwork(2).Dimensions         = 4;
  block.Dwork(2).DatatypeID         = 0;
  block.Dwork(2).Complexity         = 'Real';
  block.Dwork(2).UsedAsDiscState    = true;
  
  %stored others urc-soc
  block.Dwork(3).Name               = 'u'; 
  block.Dwork(3).Dimensions         = 2;
  block.Dwork(3).DatatypeID         = 0;
  block.Dwork(3).Complexity         = 'Real';
  block.Dwork(3).UsedAsDiscState    = true;
  
  
%endfunction

function InitConditions(block)

  %% Initialize Dwork
  %set x_prediction_k(0) and P_prediction_k(0) while calculating
  %x_k and P_k for the first iteration 
  x_p_k_0       = block.DialogPrm(2).data;
  P_p_k_0       = block.DialogPrm(3).data;

  block.Dwork(1).Data   = x_p_k_0;
  block.Dwork(2).Data   = reshape(P_p_k_0,[],1);
  block.Dwork(3).Data   = block.DialogPrm(9).data(4:5);
  
%endfunction

function Output(block)
  % assemble necessary elements
  Z_k       = block.DialogPrm(7).data;
  y_k       = block.InputPort(1).data;
  
  % state - ekf  
  x_p_k     = block.Dwork(1).Data;
  P_p_k     = reshape(block.Dwork(2).Data, 2,2);
  
  % calculate kalman-gain and real covariance matrix
  y_hat = x_p_k(2);
  C = [0, 1];
  K = P_p_k*C'/(Z_k + C*P_p_k*C');
  x_k = x_p_k + K*(y_k - y_hat);
  P_k = P_p_k - K*C*P_p_k;
  
  % assign result to Dwork
  block.Dwork(1).Data = x_k;
  block.Dwork(2).Data = reshape(P_k,[],1);
    
  % output estimated states  
  block.OutputPort(1).Data = [block.Dwork(1).Data(1)];
  
%endfunction

function Update(block)
  % update the states
  urc       = block.Dwork(3).data(2);
  soc       = block.Dwork(3).data(1);
  xk        = block.Dwork(1).Data;
  Pk        = reshape(block.Dwork(2).Data,2,2);
  Delta     = block.DialogPrm(1).Data;
  
  param     = block.DialogPrm(9).data;
  a         = block.DialogPrm(4).data;
  W_k       = block.DialogPrm(6).Data;
  Qe        = block.DialogPrm(7).Data;
  Cg        = block.DialogPrm(8).Data;
  
  % calculate prediction
  Ak = [1, 0; Delta*[param(3)/Cg,1/Delta-param(3)*((2*a(1)*soc+a(2))/Qe+...
      param(1)+1/Cg)]];
  
  x_p_kp1 = [xk(1);xk(2)+Delta*param(3)*(param(1)*param(2)*urc+xk(1)/Cg-...
      xk(2))];
  
  P_p_kp1    = Ak*Pk*Ak' + W_k;
  
  %update Dwork
  block.Dwork(1).Data   = x_p_kp1;
  block.Dwork(2).Data   = reshape(P_p_kp1,[],1);
  block.Dwork(3).Data   = [soc+Delta*xk(2)/Qe; urc+Delta*param(1)*(xk(2)-urc*param(2))];

  
%endfunction