% function [x, y, ym, L, a, t] = generate_growing_tubule(L, Nseg, delta_t, Nt, kappa, sigma_m)
%
%   Generate observations for a microtubule
%       
%       Inputs  -   L       -> initial and final lengths of microtubule
%                   lseg    -> length of each segment
%                   delta_t -> time step
%                   Nt      -> number of time steps
%                   kappa   -> bending moment for tubule
%                   gamma   -> drag coefficient
%                   sigma_m -> measurement error
%             
%       Outputs -   x -> positions down tube
%                   y(Nseg, Nt) -> displacement for each segment and
%                   timepoint
%                   ym(Nseg, Nt) -> measured displacement
%                   a(:,Nt) ->    Amplitude of each component as function
%                   of time
%                   t(Nt) -> time point

function [x, y, ym, L, a, t] = generate_growing_tubule(L, lseg, delta_t, Nt, kappa, sigma_m)

    if ( nargin < 1 )
        L = [5, 10] * 1e-6;   % Length in m
    end
    
    if (nargin < 2)
        lseg = 1.08e-7; % length of each segment
    end

    if (nargin < 3)
        delta_t = 1;  % Time steps in seconds     
    end
        
    if (nargin < 4)
        Nt = 100;  % Number of time steps
    end
 
    if (Nt>length(L))
        % Adjust length for each timestep
        minL = L(1); maxL = L(2);
        L = minL + ([1:Nt]'-1)/(Nt-1) * (maxL - minL);
    end
    
    t = [0:(Nt-1)]'*delta_t;  % Time values
        
    if  (nargin < 5)
        kappa = 2e-23  ; % Bending modulus in J.m
    end
    
    if (nargin < 6)
        sigma_m = 15e-9; % Measurement error
    end
    
    if (nargin < 7)
        gamma = 3e-3 ; % Drag coefficient (Pa.s) 
                       % Force per unit length per m/second
                       % Roughly the solution viscosity) 
    end
    
    
    Nmode = 40; % Number of modes
    a = zeros(Nmode, Nt); % Mode amplitudes at each timestep
    
    % Fill return data
    for j=1:Nt
        
        Nseg = floor(L(j)/lseg);
        x{j} = [0.5:Nseg]' * lseg;
        [Wt, qt] = Wmode(x{j},[1:Nmode],L(j), lseg, kappa); % Evaluate modes.
        tau = (gamma/kappa) * qt.^(-4); % Correlation time for each mode
        mu = exp(-delta_t./tau);     % Decay of correlation between timesteps
        dgamma = sqrt( 1- mu.*mu);    % Size of random step between each timestep.
        
        if (j==1)
           a(:,1) = randn(Nmode,1);% Initial timestep;
        else
           a(:,j) = a(:,j-1).*mu + randn(Nmode,1) .* dgamma;
        end
        
        % Calculate y
        y{j} = Wt * a(:,j);
        ym{j} = y{j} + sigma_m * randn(Nseg, 1);
    end
    
    
    

    
    
    
  