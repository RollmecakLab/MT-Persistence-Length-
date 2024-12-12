%function [y, q] = Wmode(x,n,L, delta_x, kappa) 
%   Basis of normal modes for one end fixed
%   Note modes are normalized so 
%           Integral (dy/dx)^2.dx for length = 1.
%
%   Inputs : 
%       x -> vector with center of each segment along object
%       n -> vector of mode numbers (e.g. [1 2 3] for 3 mods)
%       L -> Length of object
%       delta_x -> y(j) is the average value of mode from 
%                x(j)-delta_x/2  < x < x(j) + delta_x/2
%           If not specified, delta_x =0 and returns value at x(j)
%       kappa -> bending modulus
%
%   Outputs : 
%       y -> Array describing each mode y(j,m) is value of 
%               mode m at position x(j)
%       q -> Wave vector of each mode.
                
function [y, q] = Wmode(x,n,L, delta_x, kappa) 

    % Set segment width
    if (nargin < 4)
        delta_x = 0;
    end
    
    if (nargin < 5)
        kappa = 1; % Bending modulus
    end
    
    kT = 4e-21; % kT in J
    
    % Initialize Output
    y = zeros(length(x), length(n));
    q = zeros(length(n),1);
    
    % Process each mode.
    for j=1:length(n)

        % Normalized wave vector.
        switch n(j)
            case 1
                qn = 1.875;
            case 2 
                qn = 4.695;
            otherwise
                qn = (2*n(j) - 1) * pi/2;
        end
        
        % Wave vector
        q(j) = qn /L;
        u = x * q(j);
        du = delta_x * q(j)/2;
        
       % Version without averaging or normalization
       % y(:,j) = (cosh(qn) + cos(qn))/(sinh(qn) + sin(qn)) * (sinh(u)-sin(u)) + cos(u) - cosh(u);
       
       % Version without averaging but with normalization
       %y(:,j) = ((cosh(qn) + cos(qn))/(sinh(qn) + sin(qn)) * (sinh(u)-sin(u)) + cos(u) - cosh(u))* 2 * (L^1.5) /(qn^2);
            
       % Version with averaging and normalization   
       sf2 = ( sin(qn) - cos(qn) -  exp(-qn) ) / (1 + 2* exp(-qn) * sin(qn) - exp(-2*qn));
       sinc_val = (sin(du) + (du==0)) / (du + (du==0));
       sinhc_val = (sinh(du) + (du==0))/ (du + (du==0));
       y(:,j) =   ((cos(u) - (1-sf2 * 2 * exp(-qn)) * sin(u)) * sinc_val - ((1-sf2*exp(-qn))*exp(-u) + sf2*exp(u-qn)) * sinhc_val  ) /((L * kappa/kT)^0.5 * q(j)^2) ;

    end

