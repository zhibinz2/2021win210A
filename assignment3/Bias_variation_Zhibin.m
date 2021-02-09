%%% Assignment3 Zhibin

% Variation in beta
clear;clc;

%fix step size
dt = 0.001;
%parameter 
ndt = 0.5; %NON-DECISION TIME
mu = 0.5;  %BY CONVENTION MU MUST BE POSITIVE 
sd = 1;   %THIS IS VARIABILITY WITHIN THE WALK.  KEEP FIXED AT 1. 
nsteps = 2500; %MAX LENGTH OF WALK.  INCREASE TILL WARNING GOES AWAY
ntrials = 500; %NUMBER OF RUNS
criterion=1;


%set random number seed 
rng(19680104);


% Bias=rand(1,ntrials); % random distribution
% Bias=0.5+0.5*(0.5-Bias); % Set Mu between 0.25 to 0.75
% hist(Bias); % examine the distribution of Mu

% Bias=normrnd(0.5*dt,2*sqrt(dt),[1, ntrials]);% normal distribution
% Bias=0.5+Bias;
% hist(Bias); % examine the distribution of Mu

Bias=linspace(0.1,0.9,ntrials); % uniform distribution
% hist(Bias); % examine the distribution of Mu

% close all;

%OUTPUT VARIABLES
sample = zeros(1,nsteps+1);   %This is a single random draw from normal distribution
path = zeros(ntrials,nsteps+1); %This is all the random walks
rt = zeros(ntrials,1);  %These are the rts across trials 
correct = zeros(ntrials,1); %This is accuracy data. ZERO IS WRONG, ONE IS RIGHT


%LOOP OVER ntrials.  
for j = 1:ntrials
    goodpath = 0;
    while goodpath == 0 
        draw = normrnd(mu*dt,sd*sqrt(dt),[1,nsteps]);  %DRAW A WALK
        sample(1) = Bias(j); %START AT BIAS
        sample(2:nsteps+1) = draw; 
        walk = cumsum(sample); %SUM THE WALK.   
        crossbnd = find((walk > 1) |(walk < 0)); %TEST BOTH BOUNDARIES  
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

%compute accuracy
accuracy = mean(correct);  %COMPUTER FRACTION CORRECT

% RT for correct and error response
errorrt = rt(find(correct == 0));  % THIS IS JUST THE INCORRECT TRIALS
correctrt = rt(find(correct == 1)); %THIS IS JUST THE CORRECT TRIALS
% mean RT for correct and error response
correctrtM=median(correctrt);
errorrtM=median(errorrt); % should be fast errors


%plot all the random walks. 
subplot(1,2,2);
pathCorret=path(find(correct == 1),:);
pathError=path(find(correct == 0),:);
C=(correctrt-correctrtM).^2;
E=(errorrt-errorrtM).^2;
Pc=pathCorret(find(C==min(C)),:);
Pe=pathError(find(E==min(E)),:);


x_ax    = 0:nsteps;
X_plot  = [x_ax, fliplr(x_ax)];
Y_plot_C  = [Pc(1,:)-1.96.*sd*sqrt(dt), fliplr(Pc(1,:)+1.96.*sd*sqrt(dt))];
Y_plot_E  = [Pe(1,:)-1.96.*sd*sqrt(dt), fliplr(Pe(1,:)+1.96.*sd*sqrt(dt))];


plot(Pc(1,:)','g');
hold on;
plot(Pe(1,:)','r');
xlabel('Time (milliseconds)');
ylabel('Evidence');


fill(X_plot, Y_plot_C , 1,....
        'facecolor','green', ...
        'edgecolor','none', ...
        'facealpha', 0.3);
fill(X_plot, Y_plot_E , 1,....
        'facecolor','red', ...
        'edgecolor','none', ...
        'facealpha', 0.3);


set(gca,'YLim',[-0.5 criterion+0.5]);
set(gca,'XLim',[0 500]);
xticks(linspace(0,500,11));
xticklabels(linspace(500,1000,11));
correctTXT=['median correct RT = ' num2str(correctrtM) ' ms'];
errorTXT=['median error RT = ' num2str(errorrtM) ' ms'];

legend(correctTXT,errorTXT,'95% CI in correct path','95% CI in error path');
title('(B) median path with bias variability - FAST errors')
%Make a histogram of mean random walks. 


