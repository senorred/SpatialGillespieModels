%% Plots for protein/mRNA traces
colors=[1,0,0;0.875000000000000,1,0;0,1,0.250000000000000;0,0.625000000000000,1;0,0.250000000000000,1;1,0,1];
figure 

subplot(1,3,1)
hold on
% for k= 1:A(2)
%     plot(tvals, allprot(:,k), 'color',c(k,:))
% end
for k= 1:sets
    plot(tvals,allpgen(:,k),'color',colors(k,:),'LineWidth',3)
end
title('Protein traces','FontSize',15)
xlabel('Time','FontSize',15)
ylabel('RFU','FontSize',15)
axis([0 400 0 1600])
hold off


%% Plots for protein/mRNA noise

subplot(1,3,2)
hold on
for k= 1:A(2)
    plot(tvals, allpnoise(:,k),'color',c(k,:))
end
title('Protein Noise','FontSize',15)
ylabel('Noise','FontSize',15)
xlabel('Time','FontSize',15)
hold off


%% Plots for protein/mRNA in "noise space", CV2 vs abundance
%color modifier
    
subplot(1,3,3)
hold on
for count=1:sets
    for k= 1:trials
        %plot(mean(allpSS(k,count)),mean(allpavgcv2(k,count)),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 25);
        plot(allpSS(k,count),allpavgcv2(k,count),'color',colors(count,:),'Marker', '.','Markersize', 15);
    end
    plot(mean(allpSS(:,count)),mean(allpavgcv2(:,count)),'color',colors(count,:),'Marker', '.','Markersize', 40);
end




title('Protein CV2','FontSize',15)
set(gca,'fontsize',10);
%set(gca,'YScale','log');
%set(gca,'XScale','log');
xlabel('SS Abundance','FontSize',15)
axis([200 1800 0 0.01])
%ax = gca;
%ax.XTick = 2:2:70;
%grid on,
ylabel('CV2','FontSize',15)

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
xlabel('Time','FontSize',15)
%axis([0 400 -300 800])

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
%axis([0 400 -15000 30000])
hold off

%% Fano Factor Plots
figure
%fanovec=[0 5 10 15 20 0 20 5 15 10 0 15 5 10 20 30 25 30];
%fanovec=[0 5 10 15 20 0 20 5 15 10 0 15 5 10 20 30 25 30 0 5 10 15 20 0];
fanovec=[0 5 10 15 20 0 20 5 15 10 0 15 5 10 20 30 25 30 15];

subplot(2,1,2)
hold on
for count=1:sets
    for k= 1:trials
        plot(fanovec(count),allpfano(k,count),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 15);
        %plot(allpSS(k,count),allpfano(k,count),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 15);
    end
end
refline(0,1)
title('Protein Fano Factor','FontSize',15)
ylabel('Fano Factor','FontSize',15)
xlabel('Crowding Fraction','FontSize',15)
axis([0 30 -0.5 2.5])
%axis([200 1600 -0.5 2.5])
hold off

subplot(2,1,1)
hold on
for count=1:sets
    for k= 1:trials
        plot(fanovec(count),allmfano(k,count),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 15);
        %plot(allmSS(k,count),allmfano(k,count),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 15);
    end
end
refline(0,1)
title('mRNA Fano Factor','FontSize',15)
ylabel('Fano Factor','FontSize',15)
xlabel('Crowding Fraction','FontSize',15)
axis([0 30 -0.5 2.5])
%axis([8000 16000 -0.5 2.5])
hold off


%% CV2 vs crowding fraction
figure
%fanovec=[0 5 10 15 20 0 20 5 15 10];
%fanovec=[0 0 0 5 5 5 10 10 10 15 15 15 20 20 20];

subplot(2,1,2)
hold on
for count=1:sets
    for k= 1:trials
        plot(fanovec(count),allpavgcv2(k,count),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 15);
        %plot(allpSS(k,count),allpfano(k,count),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 15);
    end
end
refline(0,1)
title('Protein CV2','FontSize',15)
ylabel('CV2','FontSize',15)
xlabel('Crowding Fraction','FontSize',15)
set(gca,'YScale','log');
axis([0 30 0 0.01])
%axis([200 1600 -0.5 2.5])
hold off

subplot(2,1,1)
hold on
for count=1:sets
    for k= 1:trials
        plot(fanovec(count),allmavgcv2(k,count),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 15);
        %plot(allmSS(k,count),allmfano(k,count),'color',c((index(count)+(k-1)),:),'Marker', '.','Markersize', 15);
    end
end
refline(0,1)
title('mRNA CV2','FontSize',15)
ylabel('CV2','FontSize',15)
xlabel('Crowding Fraction','FontSize',15)
set(gca,'YScale','log');
axis([0 30 0 0.001])
%axis([8000 16000 -0.5 2.5])
hold off