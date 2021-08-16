function battery_dual_ekf(block)
% Level-2 MATLAB file S-Function 
% Discrete function for parameter estimation using
% extended kalman filter 
% input: battery current, output: voltage over battery and Capacitor of B-BC 
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 12; %Delta, x_0, P_0, x_L_0, P_L_0, soc-uoc_polynome, 
                        %Z_k, W_k, Z_L_k, W_L_k, Qe, Cg
  %% Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  %y(k) battery currrent
  block.InputPort(1).Dimensions        = 1;
  block.InputPort(1).DirectFeedthrough = true;
  block.InputPort(1).SamplingMode='Sample'; 
  
  %x(k)
  block.OutputPort(1).Dimensions       = [6,1];
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
  block.NumDworks = 4;
  
  %states x
  block.Dwork(1).Name               = 'x'; 
  block.Dwork(1).Dimensions         = 6;
  block.Dwork(1).DatatypeID         = 0;
  block.Dwork(1).Complexity         = 'Real';
  block.Dwork(1).UsedAsDiscState    = true;
  
  %covariance matrix P
  block.Dwork(2).Name               = 'P';
  block.Dwork(2).Dimensions         = 36;
  block.Dwork(2).DatatypeID         = 0;
  block.Dwork(2).Complexity         = 'Real';
  block.Dwork(2).UsedAsDiscState    = true;
  
  %states for iL - EKF
  block.Dwork(3).Name               = 'xL'; 
  block.Dwork(3).Dimensions         = 2;
  block.Dwork(3).DatatypeID         = 0;
  block.Dwork(3).Complexity         = 'Real';
  block.Dwork(3).UsedAsDiscState    = true;

  %states for PL - EKF
  block.Dwork(4).Name               = 'PL'; 
  block.Dwork(4).Dimensions         = 4;
  block.Dwork(4).DatatypeID         = 0;
  block.Dwork(4).Complexity         = 'Real';
  block.Dwork(4).UsedAsDiscState    = true;
%endfunction

function InitConditions(block)

  %% Initialize Dwork
  %set x_prediction_k(0) and P_prediction_k(0) while calculating
  %x_k and P_k for the first iteration 
  x_p_k_0       = block.DialogPrm(2).data;
  P_p_k_0       = block.DialogPrm(3).data;
  x_L_p_k_0     = block.DialogPrm(4).data;
  P_L_p_k_0     = block.DialogPrm(5).data;

  block.Dwork(1).Data   = x_p_k_0;
  block.Dwork(2).Data   = reshape(P_p_k_0,[],1);
  block.Dwork(3).Data   = x_L_p_k_0;
  block.Dwork(4).Data   = reshape(P_L_p_k_0,[],1);
  
  
  
%endfunction

function Output(block)
  % assemble necessary elements
  a         = block.DialogPrm(6).data;
  Z_k       = block.DialogPrm(7).data;
  y_k       = block.InputPort(1).data;
  Z_L_k     = block.DialogPrm(9).data;
  
  % state - ekf  
  x_p_k     = block.Dwork(1).Data;
  P_p_k     = reshape(block.Dwork(2).Data, 6,6);
  
  % calculate kalman-gain and real covariance matrix
  y_hat = x_p_k(6);
  C = [0, 0, 0, 0, 0, 1];
  K = P_p_k*C'/(Z_k + C*P_p_k*C');
  x_k = x_p_k + K*(y_k - y_hat);
  P_k = P_p_k - K*C*P_p_k;
  
  % i_L - ekf
  x_L_p_k   = block.Dwork(3).Data;
  P_L_p_k   = reshape(block.Dwork(4).Data, 2,2);
  
  % calculate kalman-gain and real covariance matrix
  y_L_hat = x_p_k(2);
  C_L = [0, 1];
  K_L = P_L_p_k*C_L'/(Z_L_k + C_L*P_L_p_k*C_L');
  x_L_k = x_L_p_k + K_L*(y_k-y_L_hat);
  P_L_k = P_L_p_k - K_L*C_L*P_L_p_k;
  
  % assign result to Dwork
  block.Dwork(1).Data = x_k;
  block.Dwork(2).Data = reshape(P_k,[],1);
  block.Dwork(3).Data = x_L_k;
  block.Dwork(4).Data = reshape(P_L_k,[],1);
    
  % output estimated states  
  block.OutputPort(1).Data = [block.Dwork(1).Data(1:5); block.Dwork(3).Data(1)];
  
%endfunction

function Update(block)
  % update the states
  xk        = block.Dwork(1).Data;
  Pk        = reshape(block.Dwork(2).Data,6,6);
  Delta     = block.DialogPrm(1).Data;
  xLk       = block.Dwork(3).Data;
  PLk       = reshape(block.Dwork(4).Data,2,2);
  
  a         = block.DialogPrm(6).data;
  W_k       = block.DialogPrm(8).Data;
  W_L_k     = block.DialogPrm(10).Data;
  Q_e       = block.DialogPrm(11).Data;
  Cg        = block.DialogPrm(12).Data;
  u_k       = block.Dwork(3).Data(1);
  uL_k      = block.Dwork(1).Data(5);
  
  
  % calculate prediction
  Ak = [[eye(3), zeros(3,3)]; 0, 0, 1, 0, 0, Delta/Q_e;
    Delta*(xk(6)-xk(2)*xk(5)), -Delta*xk(1)*xk(5), 0, 0, 1-Delta*xk(1)*xk(2), Delta*xk(1);
    -Delta*(xk(3)*xk(6)-xk(2)*xk(3)*xk(5)), Delta*xk(1)*xk(3)*xk(5), ...
    Delta*(u_k/Cg-xk(6)*(xk(1)+1/Cg)+xk(1)*xk(2)*xk(5)), 2*Delta*a(1)*xk(6)/Q_e,...
    Delta*xk(1)*xk(2)*xk(3), Delta*((a(2)+2*a(1)*xk(4))/Q_e-xk(3)*(xk(1)+1/Cg))+1];
  
  x_p_kp1 = [xk(1);xk(2);xk(3);xk(4)+Delta*xk(6)/Q_e;
    xk(5)+Delta*xk(1)*(xk(6)-xk(2)*xk(5)); xk(6)+Delta*(u_k*xk(3)/Cg+xk(1)*xk(2)*xk(3)*xk(5)+...
    xk(6)*((a(2)+a(1)*xk(4))/Q_e-xk(3)*(1/Cg+xk(1))))];
  
  P_p_kp1    = Ak*Pk*Ak' + W_k;
  
  ALk = [1, 0; Delta*xk(3)/Cg, 1+Delta*((a(2)+2*a(1)*xk(4))/Q_e-xk(3)*(xk(1)+1/Cg))];
  x_L_p_kp1 = [0; xLk(2)+Delta*(xLk(1)*xk(3)/Cg+xk(1)*xk(2)*xk(3)*uL_k+...
    xLk(2)*((a(2)+a(1)*xk(4))/Q_e-xk(3)*(1/Cg+xk(1))))];
  P_L_p_kp1 = ALk*PLk*ALk' + W_L_k;
 
  %update Dwork
  block.Dwork(1).Data   = x_p_kp1;
  block.Dwork(2).Data   = reshape(P_p_kp1,[],1);
  block.Dwork(3).Data   = x_L_p_kp1;
  block.Dwork(4).Data   = reshape(P_L_p_kp1,[],1);
 
%endfunction

