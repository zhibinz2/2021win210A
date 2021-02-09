%% Assignment3 Zhibin

% Variation in drift rate mu
clear;clc;

%fix step size
dt = 0.001;
%parameter 
ndt = 0.5; %NON-DECISION TIME
mu = 0.5;  %BY CONVENTION MU MUST BE POSITIVE 
sd = 1;   %THIS IS VARIABILITY WITHIN THE WALK.  KEEP FIXED AT 1. 
nsteps = 2500; %MAX LENGTH OF WALK.  INCREASE TILL WARNING GOES AWAY
ntrials = 500; %NUMBER OF RUNS
criterion = 1; %CORRECT BOUNDARY LOCATION, INCORRECT IS ZERO 
beta = 0.5; % NORMALIZED BIAS
bias = beta*criterion; %ACTUAL BIAS
%set random number seed 
rng(19680104);

nvariations=100; %BUMBER OF VARIATIONS TO DRAW FROM A DISTRIBUTION
Mu=rand(1,nvariations); % This is a vector of random draws from a random distribution of the drift rates
Accuracy=zeros(1,nvariations); % This is a vector of the acuuracy for each random drift rate
Prerror=ones(1,nvariations); % This is a vector of the probability for error response

correctrtM=zeros(1,nvariations); 
errorrtM=zeros(1,nvariations);

%LOOP OVER nvariationts of drift rate mu and then ntrials.  
for v = 1:nvariations
    
%OUTPUT VARIABLES
sample = zeros(1,nsteps+1);   %This is a single random draw from normal distribution
path = zeros(ntrials,nsteps+1); %This is all the random walks
rt = zeros(ntrials,1);  %These are the rts across trials 
correct = zeros(ntrials,1); %This is accuracy data. ZERO IS WRONG, ONE IS RIGHT  
% clear output in previous loop
clear errorrt;
clear correctrt;
    
    for j = 1:ntrials
        goodpath = 0;
        while goodpath == 0 
            draw = normrnd(Mu(v)*dt,sd*sqrt(dt),[1,nsteps]);  %DRAW A WALK
            sample(1) = bias; %START AT BIAS
            sample(2:nsteps+1) = draw; 
            walk = cumsum(sample); %SUM THE WALK.   
            crossbnd = find((walk > criterion) |(walk < 0)); %TEST BOTH BOUNDARIES  
            if ~isempty(crossbnd) %TEST IF IT CROSSED ONE OF THE BOUNDARIES AT LEAST
                goodpath = 1; %WALK IS GOOD, SET TO 1 TO EXIT WHILE LOOP
                path(j,:) = walk; %SAVE THE WALK
            else
                display('Bad Walk') %NOTIFY BAD WALK AND DRAW AGAIN LOWER
            end;
        end;
        rt(j) = crossbnd(1);  %RT IS FIRST CROSSING
        if path(j,rt(j)) > criterion  %TEST IF CORRECT
            path(j,rt(j):end) = criterion; %SET THE REST OF WALK TO BOUNDARY
            correct(j) = 1; %INDICATE CORRECT TRIAL
        else %TRIALIS INCORRECT
            path(j,rt(j):end) = 0; %SET THE REST OF WALK TO ZERO. 
        end; 
        %Add Non-decision time
        rt(j) = rt(j) + ndt/dt;	
    end
   
% compute rt mean
rtM(v)=mean(rt);


% RT for correct and error response
errorrt = rt(find(correct == 0));  % THIS IS JUST THE INCORRECT TRIALS
correctrt = rt(find(correct == 1)); %THIS IS JUST THE CORRECT TRIALS
% mean RT for correct and error response
correctrtM(v)=mean(correctrt);
errorrtM(v)=mean(errorrt);

% compute accuracy (Pr for correct)
Accuracy(v) = mean(correct);  %COMPUTER FRACTION CORRECT
% Pr for error
Prerror(v)=1-Accuracy(v);

end

% Calculate the weighted RT for correct and error response
correctrtW=zeros(1,nvariations); 
errorrtW=zeros(1,nvariations); 
% Sum the probabilities for correct and error response
sumC=sum(Accuracy);
sumE=sum(Prerror);
for v = 1:nvariations
correctrtW(v)=correctrtM(v)*Accuracy(v)/sumC;
errorrtW(v)=errorrtM(v)*Prerror(v)/sumE;
end
% Calculate Weighted mean RT
correctrtWM=sum(correctrtW); 
errorrtWM=sum(errorrtW); % should be a slow error
% correctrtWM = 759.5039
% errorrtWM = 762.2095

%%
% plot the Acuurancy against each random drift rate
figure;
subplot(2,2,1);
plot(Mu, Accuracy,'gx');
xlabel('drift rate'); ylabel('Probability');legend();
title('(A). Probability vs. drift rate variability');
hold on;
plot(Mu, Prerror, 'rx');
legend('Probability for correct response','Probability for error response');
hold off;

% plot the RT against each random drift rate
subplot(2,2,2);
plot(Mu, rtM,'x');
xlabel('drift rate'); ylabel('RT mean');
title('(B). RT mean vs. drift rate variability - slow errors');
hold on;
yline(correctrtWM,'g-','Weighted mean RT - correct','LabelVerticalAlignment','bottom');
yline(errorrtWM,'r-','Weighted mean RT - error');
label('Weighted mean RT - error')
hold off;


