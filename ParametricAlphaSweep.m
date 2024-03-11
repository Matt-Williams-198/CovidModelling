Constants;
StoredDataS = zeros(20,size(Splot,1),size(Splot,2));
StoredDataI = zeros(20,size(Splot,1),size(Splot,2));
StoredDataR = zeros(20,size(Splot,1),size(Splot,2));
count = 1;
for i = 1:20
    Di = i* 0.005;
    Covid;
    PPlotSplitter
    StoredDataS(i,:,:) = Splot(:,:).';
    StoredDataI(i,:,:) = Iplot(:,:).';
    StoredDataR(i,:,:) = Rplot(:,:).';
end
