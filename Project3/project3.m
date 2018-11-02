% Project Three : Group Three
clc;
clear;
warning off;

% Function Invoking :
enc_dec;

function enc_dec
% Loading the data file
load('data_cn_project_iii_a17.mat');

%% QUESTION 1: AUTOCORRELATION FUNCTION 
x=Stimulus;     % To store the time series
y=autocorr(x,50);   % Generate autocorrelation for tau = 0 to +50
a=autocorr(fliplr(x),50);   % Generate autocorrelation for tau = 0 to -50
z=cat(2,fliplr(a(1:50)),y);
figure;
plot([-50:50],z);     % Plotting the function
xlabel("\tau (ms)");
ylabel("R(\tau)");
title("Autocorrelation fuction of Stimulus");

%% QUESTION 2: PSTH
PSTH = zeros(4, 20000);
for i=1:4
    for j=1:50
        PSTH(i,:) = PSTH(i,:) + histcounts(All_Spike_Times{i,j}*1000,0:20000)*1000/50;
    end
end
figure;
%PSTH for 1st neuron
p1=subplot(4,1,1);
plot(PSTH(1,1:20000));
xlabel(p1,'time (ms)');
ylabel(p1,'rate (spikes/sec)');
title(['PSTH for first neuron']);
%PSTH for 2nd neuron
p2=subplot(4,1,2);
plot(PSTH(2,1:20000));
xlabel(p2,'time (ms)');
ylabel(p2,'rate (spikes/sec)');
title(['PSTH for second neuron']);
%PSTH for 3rd neuron
p3=subplot(4,1,3);
plot(PSTH(3,1:20000));
xlabel(p3,'time (ms)');
ylabel(p3,'rate (spikes/sec)');
title(['PSTH for third neuron']);
%PSTH for 4th neuron
p4=subplot(4,1,4);
plot(PSTH(4,1:20000));
xlabel(p4,'time (ms)');
ylabel(p4,'rate (spikes/sec)');
title(['PSTH for fourth neuron']);

%% QUESTION 3: POISSON PROCESS
bin_sizes = [10, 20, 50, 100, 200, 500];
for i = 1:6 % For each of the different bin sizes
    bsize = bin_sizes(i);
    figure
    
    for n=1:4   % For each of the 4 neurons
        means = [];
        variances = [];
        for j=1:50  % For each repeat for each neuron
            spike_counts = histcounts(All_Spike_Times{n,j}*1000,0:bsize:20000); % Hist of bin based counts
            M = mean(spike_counts);
            V = var(spike_counts);
            means = [means [M]];    % Append mean and variance
            variances = [variances [V]];
        end
        % Plot the means and variances for each neuron for each binsize
        % across multiple repeats of the experiment.
        max_v = max([max(means), max(variances)]);
        min_v = min([min(means), min(variances)])
        straight_line = [0,max_v];
        p = subplot(2,2,n);
        scatter(variances, means);
        hold on;
        plot(straight_line, straight_line);
        hold off;
        title(['Mean vs Variance for bin size ' num2str(bsize) 'ms - Neuron ' num2str(n)]);
        xlabel('variance');
        ylabel('mean');
        xlim([min_v max_v]);
        ylim([min_v max_v]);
    end
end
end
