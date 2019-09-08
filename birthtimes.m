birthtime=[]
for i=2:(mRNACount+1)
    birthtime(i)= find(mRNATrack(i,1,:),1,'first');  
end

