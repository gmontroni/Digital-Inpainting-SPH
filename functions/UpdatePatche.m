function BestPatche = UpdatePatche(winsize,BestPatches,nNeighbors)

BestPatche = zeros(winsize^2,1,3);

for n = 1:nNeighbors
  BestPatche = BestPatche + BestPatches(:,n,:);           %Soma os Patches.
end

BestPatche = BestPatche/nNeighbors;                       %Média.

end