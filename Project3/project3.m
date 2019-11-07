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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
       %% Question 6
    R1=corrcoef(gt1,pred1); R_squared1=R1(2)^2 
    R1
    R_squared1
    
    R2=corrcoef(gt2,pred2); R_squared2=R2(2)^2 
    R2
    R_squared2
    
    R3=corrcoef(gt3,pred3); R_squared3=R3(2)^2 
    R3
    R_squared3  
    
    R4=corrcoef(gt4,pred4); R_squared4=R4(2)^2 
    R4
    R_squared4
    
    %%the R_squared values for neuron 1 and 4 are very low and we discard
    %%them
    
    
    figure()
    subplot(2,2,1)
    scatter(gt1,pred1)
    title('PSTH vs pred for neuron 1')
    
 
    subplot(2,2,2)
    scatter(gt2,pred2)
    title('PSTH vs pred for neuron 2')
    
    
    subplot(2,2,3)
    scatter(gt3,pred3)
    title('PSTH vs pred for neuron 3')
    
    
    subplot(2,2,4)
    scatter(gt4,pred4 )
    title('PSTH vs pred for neuron 4')
    
    
    A1=zeros(100,1)
    B1=zeros(100,1)
    
    initial_r_sq=R_squared2
    old_r_sq=R_squared2
    new_r_sq=old_r_sq
    pos=zeros(1,1)
    count1=0
    while (old_r_sq - new_r_sq) < 0.01
        old_r_sq=new_r_sq
        min_=1000000000
        for i = 1:100 
            if abs(h(2,i))<min_ && abs(h(2,i))~=0
                min_=abs(h(2,i))
                pos=i
            end
        end
        if pos==0
            break
        end
        h(2,pos)=0
        pos(1,1)=0
        pred2 = conv(stimulus(15001:20000), h(2, :)); pred2 = pred2(1:5000);
        for i = 1:5000
            pred2(i) = (fit{2,1}.a)/(1+exp(-fit{2,1}.b*(pred2(i)-fit{2,1}.c)));
        end
        new_r=corrcoef(gt2,pred2); new_r_sq=new_r(2)^2
        count1=count1+1
        A1(count1,1)=count1
        B1(count1,1)=new_r_sq
        
    end
    figure()
    scatter(A1,B1)
    title('Plot of prediction performance with iterations for neuron 2')
    
    max_corr_coef_2=0
    for i =1:100
        if B1(i,1)>max_corr_coef_2
            max_corr_coef_2=B1(i,1)
        end
    end
    
    
        
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A=zeros(100,1)
    B=zeros(100,1)
    initial_r_sq=R_squared3
    old_r_sq=R_squared3
    new_r_sq=old_r_sq
    pos=zeros(1,1)
    count=0
    while (old_r_sq - new_r_sq) < 0.01
        old_r_sq=new_r_sq
        min_=1000000000
        for i = 1:100 
            if abs(h(3,i))<min_ && abs(h(3,i))~=0
                min_=abs(h(3,i))
                pos=i
            end
        end
        h(3,pos)=0
        pos(1,1)=0
        pred3 = conv(stimulus(15001:20000), h(3, :)); pred3 = pred3(1:5000);
        for i = 1:5000
            pred3(i) = (fit{3,1}.a)/(1+exp(-fit{3,1}.b*(pred3(i)-fit{3,1}.c)));
        end
        new_r=corrcoef(gt3,pred3); new_r_sq=new_r(2)^2
        count=count+1
        A(count,1)=count
        B(count,1)=new_r_sq
    end
    
    figure()
    scatter(A,B)
    title('Plot of prediction performance with iterations for neuron 3')
    max_corr_coef_3=0
    for i =1:100
        if B(i,1)>max_corr_coef_3
            max_corr_coef_3=B(i,1)
        end
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
                

    figure()
    subplot(2,2,1)
    stem(h(1,:))
    title('Linear filter for Neuron 1')
    subplot(2,2,2)
    stem(h(2,:))
    title('Linear filter for Neuron 2')
    subplot(2,2,3)
    stem(h(3,:))
    title('Linear filter for Neuron 3')
    subplot(2,2,4)
    stem(h(4,:))
    title('Linear filter for Neuron 4')
    
    
    f_t4=fft(h(4,:))
    f_t3=fft(h(3,:))
    f_t2=fft(h(2,:))
    f_t1=fft(h(1,:))
    
    vect=[-50:49]
    vect=vect'
    
    figure()
    subplot(2,2,1)
    plot(vect,f_t1)
    title('FFT of filter for Neuron 1')
    
    subplot(2,2,2)
    plot(vect,f_t2)
    title('FFT of filter for Neuron 2')
    
    subplot(2,2,3)
    plot(vect,f_t3)
    title('FFT of filter for Neuron 3')
    
    subplot(2,2,4)
    plot(vect,f_t4)
    title('FFT of filter for Neuron 4')
    
    
    fprintf('maximum predicted performance for neuron2 is %d.\n',max_corr_coef_2);
    fprintf('maximum predicted performance for neuron3 is %d.\n',max_corr_coef_3);
    
%% Part 7
%weightage of time difference
q = [0 0.001 0.01 0.1 1 10 100];
%mutual information
MI = zeros(4,100,length(q));
for iter = 1:100
    % randomly chosen 8 stimuli of length 100ms
    t = randperm(19901,8);
    for n = 1:4
        spike_segments = cell(8,50);
        for rep = 1:50
            spikes = All_Spike_Times{n,rep};
            for i = 1:8
                spike_segments{i,rep} = spikes(spikes>=t(i)/1000&spikes<(t(i)+100)/1000);
            end
        end
        
        for m = 1:length(q)
            confusion_matrix = zeros(8,8);
            for i = 1:8
                for rep = 1:50
                    mean_dist = meandist(i,rep,spike_segments,q(m));
                    mean_dist(i) = mean_dist(i)*50/49;
                    [~,k] = min(mean_dist);
                    confusion_matrix(i,k) = confusion_matrix(i,k)+1;
                end
            end
            confusion_matrix = confusion_matrix/50;
            MI(n,iter,m) = MI(n,iter,m) + MutualInfo(confusion_matrix)/100;
        end
    end
end

ci90 = c(0.9,MI);
a = MI(:,1,:);
b = ci90(:,1,:);
figure()
for n = 1:4
    subplot(2,2,n);
    plot(log10(q), a(n,:));
    hold on
    plot(log10(q), a(n,:)-b(n,:), 'r--');
    plot(log10(q), a(n,:)+b(n,:), 'r--');
    [~,p] = max(a(n,:));
    p = plot(log10(q(p)), a(n,1,p), '+');
    set(p, 'linewidth', 2)
    hold off
    title(['Discrimination - Neuron ' num2str(n)])
    xlabel('log_1_0(q)')
    ylabel('Mean MI(q) (90% confidence intervals)')
end
end

function mean_dist = meandist(i,rep,spike_segments,q)
mean_dist = zeros(1,8);
for i1 = 1:8
    for rep1 = 1:50
        if (rep1 == rep && i1 == i)
            continue
        end
        mean_dist(i1) = mean_dist(i1) + VPSDM(spike_segments{i,rep},spike_segments{i1,rep1},q);
    end
end
mean_dist = mean_dist/50;
end

function ci90 = c(ci,MI)
MIbar = mean(MI,2);
MIstd = std(abs(MI),0,2);
alpha = 1 - ci;
T_multiplier = tinv(1-alpha/2, 99);
ci90 = T_multiplier*MIstd/sqrt(99);
end

function d=VPSDM(tli,tlj,q)
nspi=length(tli);
nspj=length(tlj);

if q==0
   d=abs(nspi-nspj);
   return
elseif q==Inf
   d=nspi+nspj;
   return
end

scr=zeros(nspi+1,nspj+1);
scr(:,1)=(0:nspi)';
scr(1,:)=(0:nspj);
if(nspi && nspj)
   for i=2:nspi+1
      for j=2:nspj+1
         scr(i,j)=min([scr(i-1,j)+1 scr(i,j-1)+1 scr(i-1,j-1)+q*abs(tli(i-1)-tlj(j-1))]);
      end
   end
end
d=scr(nspi+1,nspj+1);
end

function MI = MutualInfo(confusion)
MI=0;
for i=1:size(confusion,1)
    for j=1:size(confusion,2)
        if(confusion(i,j)~=0)
            MI=MI+confusion(i,j)/size(confusion,1)*log2(confusion(i,j)/sum(confusion(:,j)));          %confusion matrix has entries of p(y/x)
        end
    end
end
end

function [fitresult, gof] = createFits(x1, y1, x2, y2, x3, y3, x4, y4)

    %% Initialization.

    % Initialize arrays to store fits and goodness-of-fit.
    fitresult = cell( 4, 1 );
    gof = struct( 'sse', cell( 4, 1 ), ...
        'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

    %% Fit: 'fit 1'.
    [xData, yData] = prepareCurveData( x1, y1 );

    % Set up fittype and options.
    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.StartPoint = [0.709754432349746 0.910192467553578 0.978691004156862];

    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft, opts );

    % Plot fit with data.
    figure( 'Name', 'fit 1' );
    h = plot( fitresult{1}, xData, yData );
    title('Fit for neuron 1')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    grid on

    %% Fit: 'fit 2'.
    [xData, yData] = prepareCurveData( x2, y2 );

    % Set up fittype and options.
    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [0.515433736542118 0.72193376784407 0.955510096008974];

    % Fit model to data.
    [fitresult{2}, gof(2)] = fit( xData, yData, ft, opts );

    % Plot fit with data.
    figure( 'Name', 'fit 2' );
    h = plot( fitresult{2}, xData, yData );
    title('Fit for neuron 2')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    grid on

    %% Fit: 'fit 3'.
    [xData, yData] = prepareCurveData( x3, y3 );

    % Set up fittype and options.
    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [0.355142740203084 0.49975884517448 0.624609065987624];

    % Fit model to data.
    [fitresult{3}, gof(3)] = fit( xData, yData, ft, opts );

    % Plot fit with data.
    figure( 'Name', 'fit 3' );
    h = plot( fitresult{3}, xData, yData );
    title('Fit for neuron 3')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    grid on

    %% Fit: 'fit 4'.
    [xData, yData] = prepareCurveData( x4, y4 );

    % Set up fittype and options.
    ft = fittype( 'a/(1+exp(-b*(x-c)))', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Robust = 'LAR';
    opts.StartPoint = [0.45087821170349 0.723648199208609 0.253875154342024];

    % Fit model to data.
    [fitresult{4}, gof(4)] = fit( xData, yData, ft, opts );

    % Plot fit with data.
    figure( 'Name', 'fit 4' );
    h = plot( fitresult{4}, xData, yData );
    title('Fit for neuron 4')
    xlabel('y(t) = s(t)*h(t)')
    ylabel('PSTH ~ \lambda(t)')
    grid on
end
