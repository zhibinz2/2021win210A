%%% Based on Deneve Latham & Pouget NN 1999's methods

clear all; close all;

% This sets up the 'grid' of the patch of cortex
totalP = 20;
P = 1:totalP;
PT = repmat(P',1,totalP);
PL = repmat(P,totalP,1);


K = 74;
C =  0.5; % contrast
v = 3.7;
sT =  0.38;
sL =  0.38; % set to be equal to sigT to match the paper
mu = 0.002;
Kw =  1;
dw =0.14+(0.718-0.14)*rand(1); % Systematically varied between 0.14 and 0.718 (p. 745)
S = 5;
sn = 5; % For Gaussian noise with fixed variance

% How many trials will we run?
N = 1000;

% Give the system an orientation and spatial frequency input
% This is in the format of an i,j coordinate, then transforms into radians
% as described just after Equation 6
T = 10; T = (2*pi*T)/totalP;
L = 10; L = (2*pi*L)/totalP;

% Initialize response network
% Here we use the same transformation as just previously, then create
% matrices out of the vectors we produce
THETA = (2*pi*P')/totalP; THETAS = repmat(THETA,1,totalP);
LAMBDA = (2*pi*P)/totalP; LAMBDAS = repmat(LAMBDA,totalP,1);

% This gives us the mean firing rate of each 'neuron' at each location
% given the current input. Using the above parameters, implement Eq 6
f = K*C*exp((cos(T-THETAS)-1)/(sT.^2)+((cos(L-LAMBDAS)-1)/(sL.^2)))+v;

for n = 1:N
    % To initialize the network, we need to add random noise to those mean
    % firing rates. This version assumes proportional noise, where variance is
    % equal to the mean because this "better approximates noise that has been 
    % measured in cortex" (p. 741).
    % That is, for a particular value of theta and lambda, total input 
    % activity aij is equal to a mean of a function (defined below) plus 
    % some noise.  The noise is related to the magnitude of the signal, 
    % s.t. aij is a drawn from a Gaussian distribution with mean 
    % f(theta,lambda) and variance also f(theta,lambda).
    % f is defined above.
    % normrnd takes in normrnd(MU,SIGMA), with SIGMA = standard deviation
    % So we must take the square root of f to provide 'sigma':
    
    alpha =  normrnd(f, sqrt(f)); % Initial conditions, see equations 6 & 7
    
    % "The initial conditions, oij(t = 0), for evolution equations 1 and
    % 2 were determined by setting oij(0) to the activity of the input
    % units, aij. The aij were obtained by drawing samples from a Gaussian 
    % distribution, as described above" (p. 741).

    o(:,:,1) = alpha; % Set up for doing this in a loop. These are the initial conditions.

    % Now compute what is going to happen on the next time step
    % (This could be done much more efficiently, but for the sake of
    % transparency we will be terrible and use nested for-loops because
    % sometimes it's easier to do that than to think in multidimensional
    % matrices)
    for i = 1:totalP
        for j = 1:totalP

            % Compute weights w in a matrix, Eq 8
            w{i,j} = Kw*exp((cos(2*pi*(THETA(i)-T)/20)-1)/dw^2+(cos(2*pi*(LAMBDA(j)-L)/20)-1)/dw^2);

            % Next, compute the us which deal with excitatory modifications of 
            % signal, and represent the activity of everybody else in the
            % network except the current neuron. Eq. 1
            u(i,j,2) = w{i,j}*alpha(i,j);

        end
    end

    % And finally, do the divisive normalization for the output on the next
    % timestep, Eq 2
    o(:,:,2) = o(:,:,1).^2./(S+mu*sum(u(:,:,2),'all'));

    % Make some plots if n = 1, to recreate Figure 2
    if n == 1 
    titles = {'Input','Output'};
    figure(1)
    for t = 1:2 
        subplot(1,2,t)
        surf(o(:,:,t))
        xlabel('Orientation (degrees)')
        ylabel('Spatial frequency (cpd)')
        % To reproduce Figure 2 from the paper, we assume the origin is
        % [0,0] and the ends of each axis are at 8cpd and 180 degrees
        set(gca,'xtick',[1,totalP],'xticklabel',[0,180],'ytick',[1,totalP],'yticklabel',[0,8])
        title(titles{t})
    end
    set(gcf,'position',[300,400,1200,600])
    end

    % Get the estimate of the location of the peak: orientation and spatial
    % frequency, so that we can get the actual theta and lambda estimates at 
    % this peak location to compute the error according to definitions in the 
    % paper for this trial of the simulation.
    % You might think it's these (just getting the location of the peak)... 
    % but it's actually not! Look at Equation 10 & the end of the Methods.
    % This is in conflict, confusingly, with what Fig 2 seems to suggest...
    % output = o(:,:,2);
    % [max_value,idx]=max(output(:));
    % [peakLoc(1),peakLoc(2)]=ind2sub(size(output),idx);
    % T_hat = THETA(peakLoc(1)); % (We don't use LAMBDA here because the thing is set to be symmetric)
    % error = (T_hat - T)^2;
    % We cannot do the above because it will depend on the discrete nature of
    % the simulation - either the peak is at the same index as the input (error
    % = 0) or it is not (error = some fixed amount). So we use their way.

    % By Equation 10, we need:
    % T_hat = phase( sum(alpha * e ^ [i*theta]) ) where i is sqrt(-1)
    % Fortunately there is a Matlab function that does this for us:
    % theta = angle(z) returns the phase angle in the interval [-?,?] for each 
    % element of a complex array z. The angles in theta are such that 
    % z = abs(z).*exp(i*theta).
    % "Each trial consisted of initializing the network with a noisy input
    % function [oij(t = 0) = aij, as described above], iterating equations 1 and 2
    % and computing ?? from equation 10"
    % Note: this seems to imply that aij is used in Eq 10, but why would we
    % need to iterate equations 1 and 2 then? So we use the output after
    % divisive normalization.
    % Sum over orientation (j) and spatial frequency (k)
    z =o(:,:,1).* exp(sqrt(-1)*THETAS); 
    T_hat = angle(z);
    error(n) = sum((T_hat - THETAS).^2,'all');

end

% Average network error
error_network = 1/(N-1) * sum(error)
