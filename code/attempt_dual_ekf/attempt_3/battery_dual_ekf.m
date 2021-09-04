function battery_dual_ekf(block)
% Level-2 MATLAB file S-Function 
% Discrete function for parameter estimation using
% extended kalman filter 
% input: battery current, output: voltage over battery and Capacitor of B-BC 
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 8; %Delta, x_0, P_0, soc-uoc_polynome, 
                        %Z_k, W_k, Qe, Cg
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
  block.OutputPort(1).Dimensions       = [8,1];
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
  block.Dwork(1).Dimensions         = 7;
  block.Dwork(1).DatatypeID         = 0;
  block.Dwork(1).Complexity         = 'Real';
  block.Dwork(1).UsedAsDiscState    = true;
  
  %covariance matrix P
  block.Dwork(2).Name               = 'P';
  block.Dwork(2).Dimensions         = 49;
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
  block.Dwork(2).Data   = reshape(P_p_k_0,[],1);
  
%endfunction

function Output(block)
  % assemble necessary elements
  a         = block.DialogPrm(4).data;
  Z_k       = block.DialogPrm(7).data;
  y_k       = block.InputPort(1).data;
  
  % state - ekf  
  x_p_k     = block.Dwork(1).Data;
  P_p_k     = reshape(block.Dwork(2).Data, 7,7);
  
  % calculate kalman-gain and real covariance matrix
  y_hat = x_p_k(7);
  C = [0, 0, 0, 0, 0, 0,1];
  K = P_p_k*C'/(Z_k + C*P_p_k*C');
  x_k = x_p_k + K*(y_k - y_hat);
  P_k = P_p_k - K*C*P_p_k;
  
  u_out = x_k(5)+a(1)*x_k(4)^2+a(2)*x_k(4)+a(3)+y_k/x_k(3);
  
  % assign result to Dwork
  block.Dwork(1).Data = x_k;
  block.Dwork(2).Data = reshape(P_k,[],1);
    
  % output estimated states  
  block.OutputPort(1).Data = [block.Dwork(1).Data; u_out];
  
%endfunction

function Update(block)
  % update the states
  xk        = block.Dwork(1).Data;
  Pk        = reshape(block.Dwork(2).Data,7,7);
  Delta     = block.DialogPrm(1).Data;
  
  a         = block.DialogPrm(4).data;
  W_k       = block.DialogPrm(6).Data;
  Qe       = block.DialogPrm(7).Data;
  Cg        = block.DialogPrm(8).Data;
  
  % calculate prediction
  Ak = [[eye(3), zeros(3,4)]; 0, 0, 0, 1, 0, 0, Delta/Qe;
    Delta*[(xk(7)-xk(2)*xk(5)), -xk(1)*xk(5), 0, 0, 1/Delta-xk(1)*xk(2), 0, xk(1);
    0, 0, 0, 0, 0, 1/Delta, 0;
    -xk(3)*(xk(7)-xk(2)*xk(5)), xk(1)*xk(3)*xk(5), ...
    xk(6)/Cg-xk(7)*(xk(1)+(a(2)+2*a(1)*xk(4))/Qe+1/Cg), -2*a(1)*xk(3)*xk(7)/Qe,...
    xk(1)*xk(2)*xk(3), xk(3)/Cg, 1/Delta-xk(3)*(xk(1)+(a(2)+2*a(1)*xk(4))/Qe+1/Cg)]];
  
  x_p_kp1 = [xk(1);xk(2);xk(3);xk(4)+Delta*xk(7)/Qe;
    xk(5)+Delta*xk(1)*(xk(7)-xk(5)*xk(2));xk(6);
    xk(7)+Delta*xk(3)*(xk(6)/Cg-xk(7)*(xk(1)+(a(2)+2*a(1)*xk(4))/Qe+1/Cg)...
    +xk(1)*xk(2)*xk(5))];
  
  P_p_kp1    = Ak*Pk*Ak' + W_k;
  
  %update Dwork
  block.Dwork(1).Data   = x_p_kp1;
  block.Dwork(2).Data   = reshape(P_p_kp1,[],1);
%endfunction

