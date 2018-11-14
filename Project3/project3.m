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
%% QUESTION 4: SPIKE TRIGGERED AVERAGE
sta = zeros(4, 100);
h = zeros(4, 100);

% Required for Whitening correction:
Rxx = autocorr(Stimulus, 99);  
Css = zeros(100,100);
Css = toeplitz(Rxx);

% Calculate the STA (<Rpx>):
for i=1:4
    total_spike_count = 0;
    for j=1:50
        % Search for how many spikes arrive before 15s (in training set)
        total_spike_count = total_spike_count + nnz(All_Spike_Times{i,j} <= 15) 
        % Generate a cummulative STA from these spikes in training set:
        for k=1:nnz(All_Spike_Times{i,j} <= 15)
            int_timing = round(All_Spike_Times{i,j}(k)*1000);
            stim_values = Stimulus(max([int_timing-99 1]):int_timing);
            partial_STA = [zeros(1,100-length(stim_values)) stim_values];
            sta(i,:) = sta(i,:) + partial_STA;  % Accumulate contribution of this repeat to sta
        end
    end
    sta(i,:)=sta(i,:)/total_spike_count; % Average out the contributions
    figure(9);
    subplot(2,2,i);
    plot(0:99,fliplr(sta(i,:)))     % Plot STA for this neuron   
    xlabel('time (ms)');
    ylabel('STA');
    ylim([-0.5 0.5]);
    title(['h(t) without correction for neuron ' num2str(i)]);
      
    % Applying whitening correction:
    h(i,100:-1:1)=(Css\sta(i,:)')';
    figure(10);
    subplot(2,2,i);
    plot(0:99,h(i,:));       % Plot corrected values
    xlabel('time (ms)');
    ylabel('h(t)');
    ylim([-1 1]);
    title(['h(t) with correction for neuron ' num2str(i)]);
end

%% Part 5
%% Q5
    stimulus = Stimulus;
    psth=PSTH/1000;
    % take first 15000 values of the stimulus to get predicted y
    pred(1, :) = conv(stimulus(1:15000), h(1, :)); x(1,1:15000) = pred(1, 1:15000); 
    pred(2, :) = conv(stimulus(1:15000), h(2, :)); x(2,1:15000) = pred(2, 1:15000);
    pred(3, :) = conv(stimulus(1:15000), h(3, :)); x(3,1:15000) = pred(3, 1:15000);
    pred(4, :) = conv(stimulus(1:15000), h(4, :)); x(4,1:15000) = pred(4, 1:15000);
    
    % estimates of the actual rates from the psth
    y(1, 1:15000) = 1000*psth(1, 1:15000); 
    y(2, 1:15000) = 1000*psth(2, 1:15000); 
    y(3, 1:15000) = 1000*psth(3, 1:15000); 
    y(4, 1:15000) = 1000*psth(4, 1:15000); 
    
    % average the results in bins
    bin_size = 30;
    for i = 1:ceil(15000/bin_size)
        e = i*bin_size;
        if e>15000
            e = 15000;
        end
        x1(i) = mean( x(1, (1+(i-1)*bin_size):e) );
        x2(i) = mean( x(2, (1+(i-1)*bin_size):e) );
        x3(i) = mean( x(3, (1+(i-1)*bin_size):e) );
        x4(i) = mean( x(4, (1+(i-1)*bin_size):e) );
    
        y1(i) = mean( y(1, (1+(i-1)*bin_size):e) );
        y2(i) = mean( y(2, (1+(i-1)*bin_size):e) );
        y3(i) = mean( y(3, (1+(i-1)*bin_size):e) );
        y4(i) = mean( y(4, (1+(i-1)*bin_size):e) );
    end
    
    figure()
    subplot(2,2,1)
    scatter(x1, y1)
    title('PSTH vs y(t) for neuron 1')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    subplot(2,2,2)
    scatter(x2, y2)
    title('PSTH vs y(t) for neuron 2')
    xlabel('y(t) = s(t)*h(t)')  
    ylabel('PSTH ~ \lambda(t)')
    subplot(2,2,3)
    scatter(x3, y3)
    title('PSTH vs y(t) for neuron 3')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    subplot(2,2,4)
    scatter(x4, y4)
    title('PSTH vs y(t) for neuron 4')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
        
    [fit, gof] = createFits(x1, y1, x2, y2, x3, y3, x4, y4);
    
    save('fits', 'fit', 'gof')
    
    % predictions
    % linear filter
    pred1 = conv(stimulus(15001:20000), h(1, :)); pred1 = pred1(1:5000);
    pred2 = conv(stimulus(15001:20000), h(2, :)); pred2 = pred2(1:5000);
    pred3 = conv(stimulus(15001:20000), h(3, :)); pred3 = pred3(1:5000);
    pred4 = conv(stimulus(15001:20000), h(4, :)); pred4 = pred4(1:5000);
    % add non-linearity
    for i = 1:5000
        pred1(i) = (fit{1,1}.a)/(1+exp(-fit{1,1}.b*(pred1(i)-fit{1,1}.c)));
        pred2(i) = (fit{2,1}.a)/(1+exp(-fit{2,1}.b*(pred2(i)-fit{2,1}.c)));
        pred3(i) = (fit{3,1}.a)/(1+exp(-fit{3,1}.b*(pred3(i)-fit{3,1}.c)));
        pred4(i) = (fit{4,1}.a)/(1+exp(-fit{4,1}.b*(pred4(i)-fit{4,1}.c)));
    end 
    
    % ground truths
    gt1 = 1000*psth(1, 15001:20000);
    gt2 = 1000*psth(2, 15001:20000);
    gt3 = 1000*psth(3, 15001:20000);
    gt4 = 1000*psth(4, 15001:20000);
end
