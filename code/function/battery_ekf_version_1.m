function battery_ekf_version_1(block)
% Level-2 MATLAB file S-Function 
% Discrete-time state-space model with u and 
% a (time-varying parameter) as inputs
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 7; %Delta, x0, P0, soc-uoc_polynome, 
                        %Z_k, W_k, Qe
  %% Register number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  %u(k)
  block.InputPort(1).Dimensions        = 1;
  block.InputPort(1).DirectFeedthrough = true;
  block.InputPort(1).SamplingMode='Sample'; 
  %y(k)
  block.InputPort(2).Dimensions        = 1;
  block.InputPort(2).DirectFeedthrough = true;
  block.InputPort(2).SamplingMode='Sample'; 
  
  %x(k)
  block.OutputPort(1).Dimensions       = [5,1];
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
  block.Dwork(1).Dimensions         = 5;
  block.Dwork(1).DatatypeID         = 0;
  block.Dwork(1).Complexity         = 'Real';
  block.Dwork(1).UsedAsDiscState    = true;
  
  %covariance matrix P
  block.Dwork(2).Name               = 'P';
  block.Dwork(2).Dimensions         = 25;
  block.Dwork(2).DatatypeID         = 0;
  block.Dwork(2).Complexity         = 'Real';
  block.Dwork(2).UsedAsDiscState    = true;

%endfunction

function InitConditions(block)

  %% Initialize Dwork
  %set x_prediction_k(0) and P_prediction_k(0) while calculating
  %x_k and P_k for the first iteration 
  x_p_k_0     = block.DialogPrm(2).data;
  P_p_k_0     = block.DialogPrm(3).data;

  block.Dwork(1).Data   = x_p_k_0;
  block.Dwork(2).Data   = reshape(P_p_k_0,[],1);
  
  
%endfunction

function Output(block)
  %assemble necessary elements
  a         = block.DialogPrm(4).data;
  Z_k       = block.DialogPrm(5).data;
  u_k       = block.InputPort(1).data;
  y_k       = block.InputPort(2).data;
  
  x_p_k     = block.Dwork(1).Data;
  P_p_k     = reshape(block.Dwork(2).Data, 5,5);
  
  %calculate kalman-gain and real covariance matrix
  y_hat = a(1)*x_p_k(4)^2+a(2)*x_p_k(4)+a(3)+x_p_k(5)+...
      u_k/x_p_k(3);
  C = [0, 0, -u_k/x_p_k(3)^2, 2*a(1)*x_p_k(4)+a(2), 1];
  K = P_p_k*C'/(Z_k + C*P_p_k*C');
  xk = x_p_k + K*(y_k - y_hat);
  Pk = P_p_k - K*C*P_p_k;
  
  %assign result to Dwork
  block.Dwork(1).Data = xk;
  block.Dwork(2).Data = reshape(Pk,[],1);
    
  %output estimated states  
  block.OutputPort(1).Data = block.Dwork(1).Data;
  
%endfunction

function Update(block)
  %update the states
  xk       = block.Dwork(1).Data;
  Pk       = reshape(block.Dwork(2).Data,5,5);
  Delta     = block.DialogPrm(1).Data;
  W_k       = block.DialogPrm(6).Data;
  Q_e       = block.DialogPrm(7).Data;
  
  u_k       = block.InputPort(1).Data;
  
  %calculate prediction
  Ak       = [[eye(4); zeros(1,4)], [Delta*(u_k-xk(5)*xk(3));0;...
      -Delta*xk(1)*xk(5);0;1-Delta*xk(1)*xk(3)]];
    
  x_p_kp1    = [xk(1);xk(2);xk(3);xk(4)+Delta*u_k/Q_e;...
      xk(5)-Delta*xk(1)*(xk(5)*xk(2)-u_k)];
  P_p_kp1    = Ak*Pk*Ak' + W_k;
 
  %update Dwork
  block.Dwork(1).Data   = x_p_kp1;
  block.Dwork(2).Data   = reshape(P_p_kp1,[],1);
 
%endfunction
