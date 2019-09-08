%% Master Plot for SXPM 2.0 trials
% Description: Wrapper for plotting multiple simulation runs
% 1000 second runs, 50 iterations per condition
% 20x20 flat XY area
% 400 Ribosomes, seeded randomly
% alpha= kCM*0.001
% dRibo= 1.0
% dPoly= 0.0001
% Reflective Boundaries

%% MODEL -- SXPM 2.0

alphaConditions= [0.001];
crowdConditions= [0, 10, 20, 30, 40, 50];
iterations= 1;
figure
hold on

for I=1:length(alphaConditions)
    II= alphaConditions(I)
    for J= 1:length(crowdConditions)
        JJ= crowdConditions(J)
%         figure
%         hold on
        for K= 1:iterations
            KK= K;
            
            c=hsv(length(crowdConditions));
            
            
%% clear variables, seed, and name the run
clearvars -except alphaConditions crowdConditions iterations II JJ KK I J K c

Name = sprintf('alpha%gCrowd%giter%g.mat',II,JJ,KK);
load(Name)


plot(tspan(1:1000), mRNATrack2(1:1000), 'color', c(J,:))
title('mRNA Count over time','FontSize',15)
set(gca,'fontsize',10);
%set(gca,'YScale','log');
%set(gca,'XScale','log');
xlabel('Time','FontSize',15)
%axis([200 1800 0 0.01])
%ax = gca;
%ax.XTick = 2:2:70;
%grid on
ylabel('mRNA Count','FontSize',15)


        end
        %hold off
    end
end
hold off
