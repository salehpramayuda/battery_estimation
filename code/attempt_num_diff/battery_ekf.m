function battery_ekf(block)
% Level-2 MATLAB file S-Function 
% Discrete function for parameter estimation using extended kalman filter 
% Input: battery current; Output: voltage over battery (and parallel
% capacitor)
  setup(block)
  
%endfunction

function setup(block)
  
% Register number of dialog parameters
  block.NumDialogPrms = 8; %Delta, x0, P0, a, 
                        %Z_k, W_k, Qe, Cg
  %% Register number of input and output ports
  block.NumInputPorts  = 3;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
  
  % u(k) battery current
  block.InputPort(1).Dimensions         = 1;
  block.InputPort(1).DirectFeedthrough  = true;
  block.InputPort(1).SamplingMode       ='Sample'; 
  
  % y(k) measured voltage
  block.InputPort(2).Dimensions         = 1;
  block.InputPort(2).DirectFeedthrough  = true;
  block.InputPort(2).SamplingMode       ='Sample'; 
  
  % y(k+1) measured voltage
  block.InputPort(3).Dimensions         = 1;
  block.InputPort(3).DirectFeedthrough  = true;
  block.InputPort(3).SamplingMode       ='Sample'; 
  
  % x(k) estimated states
  block.OutputPort(1).Dimensions        = 6;
  block.OutputPort(1).SamplingMode      ='Sample';
  
%   % y_hat(k) 'simulated' voltage
%   block.OutputPort(2).Dimensions        = 1;
%   block.OutputPort(2).SamplingMode      = 'Sample';
  
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
  
  % States x
  block.Dwork(1).Name               = 'x'; 
  block.Dwork(1).Dimensions         = 6;
  block.Dwork(1).DatatypeID         = 0;
  block.Dwork(1).Complexity         = 'Real';
  block.Dwork(1).UsedAsDiscState    = true;
  
  % Covariance matrix P
  block.Dwork(2).Name               = 'P';
  block.Dwork(2).Dimensions         = 36;
  block.Dwork(2).DatatypeID         = 0;
  block.Dwork(2).Complexity         = 'Real';
  block.Dwork(2).UsedAsDiscState    = true;

%endfunction

function InitConditions(block)

  %% Initialize Dwork
  % Set x_prediction_k(0) and P_prediction_k(0) while calculating
  % x_k and P_k for the first iteration 
  x_p_k_0       = block.DialogPrm(2).data;
  P_p_k_0       = block.DialogPrm(3).data;

  block.Dwork(1).data   = x_p_k_0;
  block.Dwork(2).data   = reshape(P_p_k_0,[],1);
%endfunction

function Output(block)
  % Assemble necessary elements
  Z_k   = block.DialogPrm(5).data;
  y_k   = block.InputPort(2).data;
  
  x_p_k = block.Dwork(1).data;
  P_p_k = reshape(block.Dwork(2).data, length(x_p_k),length(x_p_k));
  
  % Calculate kalman-gain and real covariance matrix
  y_hat = x_p_k(6); % ib*R0 + u_oc + u_rc
  C     = [0, 0, 0, 0, 0, 1];
  
  K     = P_p_k*C'/(Z_k + C*P_p_k*C');
  x_k   = x_p_k + K*(y_k - y_hat);
  P_k   = P_p_k - K*C*P_p_k;
  
  % Assign result to Dwork
  block.Dwork(1).data = x_k;
  block.Dwork(2).data = reshape(P_k,[],1);
    
  % Output estimated states  
  block.OutputPort(1).data = block.Dwork(1).data;
  
%endfunction

function Update(block)
  % Access data for state update
  Delta = block.DialogPrm(1).data;
  a     = block.DialogPrm(4).data;
  W_k   = block.DialogPrm(6).data;
  Qe    = block.DialogPrm(7).data;
  Cg    = block.DialogPrm(8).data;
  uk    = block.InputPort(1).data;
  ykp1  = block.InputPort(3).data;
  xk    = block.Dwork(1).data;
  Pk    = reshape(block.Dwork(2).data,length(xk),length(xk));
  
  i_c   = Cg*(ykp1 - uk/xk(3)-xk(5)-a(1)*xk(4)^2-a(2)*xk(4)-a(3))/Delta;
  
  % Calculate prediction for state and cov-Matrix
  Ak        = [[eye(4), zeros(4,2)]; Delta*[uk-xk(2)*xk(5), -xk(1)*xk(5),...
      0, 0, 1/Delta - xk(1)*xk(2),0]; [0,0, uk/xk(3)^2, -2*a(1)*xk(4)-a(2), -1, 1]];
  x_p_kp1   = [xk(1); xk(2); xk(3); xk(4)+Delta*uk/Qe; 
      xk(5)+Delta*xk(1)*(uk-xk(2)*xk(5)); xk(6)+Delta*i_c/Cg];
  
  P_p_kp1    = Ak*Pk*Ak' + W_k;
 
  % Update Dwork
  block.Dwork(1).data   = x_p_kp1;
  block.Dwork(2).data   = reshape(P_p_kp1,[],1);
%endfunction
