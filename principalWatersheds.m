clear
load watershedHistogram.mat

threshold= .80;
prinMatrix= zeros(iterations,length(areaConditions));

for i= 1:length(areaConditions)
    for j= 1:iterations
        AA=areaWinnersStore(:,j,i);
        BB= AA(~isnan(AA))/sum(AA,'omitnan');
        CC= cumsum(sort(BB,'descend'));
        principals= find(CC>threshold, 1);
        prinMatrix(j,i)= principals;
    end
end

%histogram(prinMatrix)
figure
hist(prinMatrix,50)