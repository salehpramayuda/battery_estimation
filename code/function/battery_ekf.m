function battery_ekf(block)
% Level-2 MATLAB file S-Function 
% Discrete function for parameter estimation using
% extended kalman filter 
% input: battery current, output: voltage over battery and Capacitor of B-BC 
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 8; %Delta, x0, P0, soc-uoc_polynome, 
                        %Z_k, W_k, Qe, u_km1
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
  block.NumDworks = 3;
  
  %states x
  block.Dwork(1).Name               = 'x'; 
  block.Dwork(1).Dimensions         = 8;
  block.Dwork(1).DatatypeID         = 0;
  block.Dwork(1).Complexity         = 'Real';
  block.Dwork(1).UsedAsDiscState    = true;
  
  %covariance matrix P
  block.Dwork(2).Name               = 'P';
  block.Dwork(2).Dimensions         = 64;
  block.Dwork(2).DatatypeID         = 0;
  block.Dwork(2).Complexity         = 'Real';
  block.Dwork(2).UsedAsDiscState    = true;
  
  %previous voltage u_km1
  block.Dwork(3).Name               = 'u_km1'; 
  block.Dwork(3).Dimensions         = 1;
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
  u_km1         = block.DialogPrm(8).data;

  block.Dwork(1).Data   = x_p_k_0;
  block.Dwork(2).Data   = reshape(P_p_k_0,[],1);
  block.Dwork(3).Data   = u_km1;
  
  
  
%endfunction

function Output(block)
  %assemble necessary elements
  a         = block.DialogPrm(4).data;
  Z_k       = block.DialogPrm(5).data;
  u_k       = block.InputPort(1).data;
  y_k       = block.InputPort(2).data;
  
  x_p_k     = block.Dwork(1).Data;
  P_p_k     = reshape(block.Dwork(2).Data, 8,8);
  
% calculate kalman-gain and real covariance matrix

  u_oc = a(1)*x_p_k(4)^2+a(2)*x_p_k(4)+a(3);
  y_hat = x_p_k(4)*x_p_k(8);
  C = [0, 0, 0, x_p_k(8), 0, 0, 0, x_p_k(4)];
  
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
  xk        = block.Dwork(1).Data;
  Pk        = reshape(block.Dwork(2).Data,8,8);
  Delta     = block.DialogPrm(1).Data;
  a         = block.DialogPrm(4).data;
  W_k       = block.DialogPrm(6).Data;
  Q_e       = block.DialogPrm(7).Data;
  u_km1     = block.Dwork(3).Data;
  u_k       = block.InputPort(1).Data;
  gamma     = u_km1/u_k;
  
  
  %calculate prediction
  ic = xk(7)-2e-2-u_k;
  Ak = [[eye(5) zeros(5,3)];...
      Delta*(u_k-xk(2)*xk(6)),-Delta*xk(1)*xk(6),0,0,0,1-Delta*xk(1)*xk(2),0,0;...
      -Delta*xk(3)*(u_k-xk(2)*xk(6)), Delta*xk(1)*xk(3)*xk(6), -Delta*(xk(4)*...
      (-ic) + xk(1)*(u_k-xk(2)*xk(6)) + u_k*(a(2)+2*a(1)*xk(5))/Q_e), -Delta*...
      (xk(3)*(-ic)-(gamma-1)*ic/Delta), -2*Delta*a(1)*u_k*xk(3)/Q_e, Delta*...
      xk(1)*xk(2)*xk(3), Delta*(xk(3)*xk(4)-xk(4)*(gamma-1)/Delta), 0;...
      0, 0, 0, 0, 0, 0, Delta, 1];
  
  x_p_kp1 = [xk(1);xk(2);xk(3); xk(4); xk(5)+Delta*u_k/Q_e; xk(6)+Delta*xk(1)*(u_k-...
      xk(2)*xk(6)); xk(7)-Delta*(xk(3)*(xk(4)*(-ic)+xk(1)*(u_k-xk(2)*xk(6))+ ...
      u_k/Q_e*(a(2)+2*a(1)*xk(5)))+xk(4)*(gamma-1)*ic/Delta); xk(8)+Delta*ic];
  
  P_p_kp1    = Ak*Pk*Ak' + W_k;
 
  %update Dwork
  block.Dwork(1).Data   = x_p_kp1;
  block.Dwork(2).Data   = reshape(P_p_kp1,[],1);
  block.Dwork(3).Data   = u_k;
 
%endfunction
