function battery_ekf(block)
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
  block.InputPort(1).DirectFeedthrough = false;
  block.InputPort(1).SamplingMode='Sample'; 
  %y(k)
  block.InputPort(2).Dimensions        = 1;
  block.InputPort(2).DirectFeedthrough = false;
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
  x_pk0     = block.DialogPrm(2).data;
  P_pk0     = block.DialogPrm(3).data;

  block.Dwork(1).Data   = x_pk0;
  block.Dwork(2).Data   = reshape(P_pk0,[],1);
  
  
%endfunction

function Output(block)
  %assemble necessary elements
  a         = block.DialogPrm(4).data;
  Z_k       = block.DialogPrm(5).data;
  u_k       = block.InputPort(1).data;
  y_k       = block.InputPort(2).data;
  
  x_pk     = block.Dwork(1).Data;
  P_pk     = reshape(block.Dwork(2).Data, 5,5);
  
  %calculate kalman-gain and real covariance matrix
  y_hat_k = a(1)*x_pk(4)^2+a(2)*x_pk(4)+a(3)+x_pk(5)+...
      u_k/x_pk(3);
  C_k = [0, 0, -u_k/x_pk(3)^2, 2*a(1)*x_pk(4)+a(2), 1];
  K_k = P_pk*C_k'/(Z_k + C_k*P_pk*C_k');
  x_k = x_pk + K_k*(y_k - y_hat_k);
  P_k = P_pk - K_k*C_k*P_pk;
  
  %assign result to Dwork
  block.Dwork(1).Data = x_k;
  block.Dwork(2).Data = reshape(P_k,[],1);
    
  %output estimated states  
  block.OutputPort(1).Data = block.Dwork(1).Data;
  
%endfunction

function Update(block)
  %update the states
  x_k       = block.Dwork(1).Data;
  P_k       = reshape(block.Dwork(2).Data,5,5);
  Delta     = block.DialogPrm(1).Data;
  W_k       = block.DialogPrm(6).Data;
  Q_e       = block.DialogPrm(7).Data;
  
  u_k       = block.InputPort(1).Data;
  
  %calculate prediction
  A_k       = [[eye(4); zeros(1,4)], [Delta*(u_k-x_k(5)*x_k(3));0;...
      -Delta*x_k(1)*x_k(5);0;1-Delta*x_k(1)*x_k(3)]];
    
  x_pkp1    = [x_k(1);x_k(2);x_k(3);x_k(4)+Delta*u_k/Q_e;...
      x_k(5)-Delta*x_k(1)*(x_k(5)*x_k(2)-u_k)];
  P_pkp1    = A_k*P_k*A_k' + W_k;
 
  %update Dwork
  block.Dwork(1).Data   = x_pkp1;
  block.Dwork(2).Data   = reshape(P_pkp1,[],1);
 
%endfunction
