%% Noise analysis tool for multi-condition protein/mRNA traces with different sizes


%% Variables & Preallocated Matrices for all sets

% mtraces= ;    % plug these in beforehand into the .mat file
% ptraces= ;    % these matrices should ALWAYS be the same size
                % if not, ya dun goofed
% labels= {'label1', 'label2', 'label3'};   % legend labels for each set

allprot= ptraces;
allmrna= mtraces;
A= size(allprot);       % A(1)= number of timepoints, A(2)= sets*trials

% Specify set number and trials within each set, if necessary. 
% Try to just save this information in the .mat file instead
% sets= ;
% trials= ;

% Define variables used for ALL sets
% Protein
allprotzero= zeros(A);
allpSS= zeros(trials, sets);        % Steady State
allpgen= zeros(A(1),sets);          % General Trends for zeroed trace
allpgengraph= zeros(A(1),sets);     % General Trends
allpnoise= zeros(A);                % Noise Traces
allpauto= zeros((A(1)*2)-1, A(2));  % Autocorrelation Traces
allpautoavg= zeros((A(1)*2)-1,sets);% Average Autocorrelation Traces
allpvar= zeros(trials, sets);       % Variance
allpavgcv2= zeros(trials, sets);    % CV2
allpfano= zeros(trials,sets);       % Fano factor

% mRNA
allmrnazero= zeros(A);
allmSS= zeros(trials, sets);        % Steady State
allmgen= zeros(A(1),sets);          % General Trends for zeroed trace
allmgengraph= zeros(A(1),sets);     % General Trends
allmnoise= zeros(A);                % Noise Traces
allmauto= zeros((A(1)*2)-1, A(2));  % Autocorrelation Traces
allmautoavg= zeros((A(1)*2)-1,sets);% Average Autocorrelation Traces
allmvar= zeros(trials, sets);       % Variance
allmavgcv2= zeros(trials, sets);    % CV2
allmfano= zeros(trials,sets);       % Fano factor

%% Loop through individual sets of the same number of traces

for count=1:sets
    count
    % Define variables within an individual set
    index=1:trials:A(2);
    ptraces=[];
    mtraces=[];
    ptraces=allprot(:,index(count):index(count)+(trials-1));    %choose first set
    mtraces=allmrna(:,index(count):index(count)+(trials-1));

    S= size(ptraces);
    numtraces=S(2);     % number of traces in each set (should be same as "trials")
    timepoints=S(1);    % number of timepoints in each trace, same as A(1)

    SSptraces= zeros(numtraces,1);          % Steady State for protein traces
    SSmtraces= zeros(numtraces,1);          % Steady State for mRNA traces
    GainStore= zeros(numtraces,1);          % Found Gain Storage for protein traces
    mGainStore= zeros(numtraces,1);         % Found Gain Storage for mRNA traces
    pnoise= zeros(timepoints,numtraces);    % Noise traces for protein
    pauto= zeros(timepoints*2-1,numtraces); % Autocorrelation traces for protein
    pvar= zeros(numtraces,1);               % Variance for protein traces
    pavgcv2= zeros(numtraces,1);            % CV2 for protein traces
    pfano= zeros(numtraces,1);
    mnoise= zeros(timepoints,numtraces);    % "" for mRNA
    mauto= zeros(timepoints*2-1,numtraces); %
    mvar= zeros(numtraces,1);               %
    mavgcv2= zeros(numtraces,1);            %
    mfano= zeros(numtraces,1);
    tvals= 0:5:(timepoints-1)*5;            % Timepoints in minutes, starting from 0
    
    % calculate protein and mrna mean trace for the non-normalized graph

    pavggraph= nanmean(ptraces,2);
    mavggraph= nanmean(mtraces,2);

    % Subtract first value from each trace
    % may later want to remove the lowest value from the mRNA traces then
    % subtract the first value and truncate the corresponding protein traces
    % This is the first step in normalization.

    for i= 1:numtraces
        ptraces(:,i)= ptraces(:,i)-ptraces(1,i);
        SSptraces(i)= ptraces(end,i);
        mtraces(:,i)= mtraces(:,i)-mtraces(1,i);
        SSmtraces(i)= mtraces(end,i);
    end

    % calculate protein and mrna mean trace

    pavg= nanmean(ptraces,2);
    mavg= nanmean(mtraces,2);

    % Coarse grain gain finder for protein and mRNA (Thanks, Chuck!)

    for j=1:numtraces

        % Generate the Noise traces
        pnoise(:,j)= ptraces(:,j)- pavg;     % protein
        mnoise(:,j)= mtraces(:,j)- mavg;   % mrna

        % Generate the Autocorrelations and Variance
        % Protein
        pautoTemp= xcorr(pnoise(:,j),'biased');
        pauto(:,j)=pautoTemp;
        pvar(j)= pautoTemp(timepoints);
        pavgcv2(j)= pautoTemp(timepoints)/(SSptraces(j))^2;
        pfano(j)= pautoTemp(timepoints)/(SSptraces(j));
        % mRNA
        mautoTemp= xcorr(mnoise(:,j),'biased');
        mauto(:,j)=mautoTemp;
        mvar(j)= mautoTemp(timepoints);
        mavgcv2(j)= mautoTemp(timepoints)/(SSmtraces(j))^2;
        mfano(j)= mautoTemp(timepoints)/(SSmtraces(j));

    end
    
    % store the steady state, noise, variance, and cv2 for this particular
    % set (each "count" iteration) in the matrices defined for all sets

    
    % Protein
    allprotzero(:,index(count):index(count)+(trials-1))= ptraces;
    allpSS(:,count)= SSptraces;
    allpgen(:,count)= pavg;
    allpgengraph(:,count)= pavggraph;
    allpnoise(:,index(count):index(count)+(trials-1))= pnoise;
    allpauto(:,index(count):index(count)+(trials-1))= pauto;
    allpautoavg(:,count)= mean(pauto,2);
    allpvar(:,count)= pvar;
    allpavgcv2(:,count)= pavgcv2;
    allpfano(:,count)=pfano;
    % mRNA
    allmrnazero(:,index(count):index(count)+(trials-1))= mtraces;
    allmSS(:,count)= SSmtraces;
    allmgen(:,count)= mavg;
    allmgengraph(:,count)= mavggraph;
    allmnoise(:,index(count):index(count)+(trials-1))= mnoise;
    allmauto(:,index(count):index(count)+(trials-1))= mauto;
    allmautoavg(:,count)= mean(mauto,2);
    allmvar(:,count)= mvar;
    allmavgcv2(:,count)= mavgcv2;
    allmfano(:,count)=mfano;
    
end

%% Color all trials in each set of conditions the same color along hsv spectrum
% This makes the traces into pretty pretty rainbows! Hooray!

colors= hsv(sets);
c=[];
for cc= 1:sets
    b=zeros(trials,3);
    for ccc=1:trials
        b(ccc,:)=colors(cc,:);
    end
    c= vertcat(c,b);    % concatenating in a nested for loop #livinontheedge
end

%c = hsv(A(2)); % this just does every trace an individual color

%% Plots for protein/mRNA traces
figure 

subplot(2,3,4)
hold on
for k= 1:A(2)
    plot(tvals, allprot(:,k), 'color',c(k,:))
end
for k= 1:sets
    plot(tvals,allpgengraph(:,k),'color',colors(k,:),'LineWidth',3)
end
title('Protein traces','FontSize',15)
xlabel('Time (min)','FontSize',15)
ylabel('RFU','FontSize',15)
hold off

subplot(2,3,1)
hold on
for k= 1:A(2)
    plot(tvals, allmrna(:,k),'color',c(k,:))
end
for k= 1:sets
    plot(tvals,allmgengraph(:,k),'color',colors(k,:),'LineWidth',3)
end
title('mRNA traces','FontSize',15)
xlabel('Time (min)','FontSize',15)
ylabel('RFU','FontSize',15)
hold off

%% Plots for protein/mRNA noise

subplot(2,3,5)
hold on
for k= 1:A(2)
    plot(tvals, allpnoise(:,k),'color',c(k,:))
end
title('Protein Noise','FontSize',15)
ylabel('Noise','FontSize',15)
xlabel('Time (min)','FontSize',15)
hold off

subplot(2,3,2)
hold on
for k= 1:A(2)
    plot(tvals, allmnoise(:,k),'color',c(k,:))
end
title('mRNA Noise','FontSize',15)
ylabel('Noise','FontSize',15)
xlabel('Time (min)','FontSize',15)
hold off

%% Plots for protein/mRNA in "noise space", CV2 vs abundance
    
subplot(2,3,6)
hold on
for count=1:sets
    for k= 1:trials
        plot(allpSS(k,count),allpavgcv2(k,count),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 15);
    end
end
title('Protein CV2','FontSize',15)
set(gca,'fontsize',10);
set(gca,'YScale','log');
set(gca,'XScale','log');
xlabel('SS Abundance','FontSize',15)
axis([1 10000 0.0001 1])
%ax = gca;
%ax.XTick = 2:2:70;
%grid on
ylabel('CV2','FontSize',15)

hold off

subplot(2,3,3)
hold on
% Get legend values
LEG=zeros(sets,1);
for count=1:sets
    leg=zeros(trials,1);
    for k= 1:trials
        leg(k)=plot(allmSS(k,count),allmavgcv2(k,count),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 15);
    end
    LEG(count)=leg(1);
end
title('mRNA CV2','FontSize',15)
legend(LEG,labels);
set(legend,'Position',[0.9 0.3 0.1 0.5]);
set(gca,'fontsize',10);
set(gca,'YScale','log');
set(gca,'XScale','log');
axis([100 100000 0.00001 0.1])
xlabel('SS Abundance','FontSize',15)
ylabel('CV2','FontSize',15)
%ax = gca;
%ax.XTick = 2:2:70;
%grid on

hold off

%% Autocorrelation plots
figure
subplot(2,1,2)
hold on
for k= 1:A(2)
    plot(tvals, allpauto(A(1):end,k),'color',c(k,:))
end
for k= 1:sets
    plot(tvals,allpautoavg(A(1):end,k),'color',colors(k,:),'LineWidth',3)
end
title('Protein Autocorrelation Traces','FontSize',15)
ylabel('Autocorrelation','FontSize',15)
xlabel('Time (min)','FontSize',15)

hold off

subplot(2,1,1)
hold on
for k= 1:A(2)
    plot(tvals, allmauto(A(1):end,k),'color',c(k,:))
end
for k= 1:sets
    plot(tvals,allmautoavg(A(1):end,k),'color',colors(k,:),'LineWidth',3)
end
title('mRNA Autocorrelation Traces','FontSize',15)
ylabel('Autocorrelation','FontSize',15)
xlabel('Time (min)','FontSize',15)
hold off

%%