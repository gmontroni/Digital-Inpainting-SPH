function [DcompBinaryList, R, Dcomp, Chosendist] = ... 
DistComp(distPatche, BinaryMask, Patches, sizeRows, sizeCols, winsize, iCoordsVSample, jCoordsVSample, halfwinsize, nNeighbors)

  [nPixelsPerPatch, nPatches, nChannels] = size(Patches);            %Dimensão do conjunto Patches
  limitThreshold = 0.3;                                              
  MaskRow = BinaryMask(:)';                                          %Vetorr coluna da máscara.
  nDist = zeros(nPixelsPerPatch,2,nPatches);                         %Vetor de coordenadas de cada Patche.
  
  %--------------------------------------------------
  % STEP 1.
  %--------------------------------------------------
 
  for n = 1 : nPatches
   
  iRange = max(iCoordsVSample(n)-halfwinsize,1):min(iCoordsVSample(n)+halfwinsize,sizeRows);       %tamanho da linha
  jRange = max(jCoordsVSample(n)-halfwinsize,1):min(jCoordsVSample(n)+halfwinsize,sizeCols);       %tamanho da coluna
  
  iAuxDist = iRange'*ones(1,winsize);
  jAuxDist = jRange'*ones(1,winsize);
  
  nDist(:,1,n) = iAuxDist(:); 
  nDist(:,2,n) = sort(jAuxDist(:));
  
  end
  
  %--------------------------------------------------
  % STEP 2.
  %--------------------------------------------------
  Dcomp = zeros(1,nPatches);
  
  for n = 1 : nPatches
    %Dcomp = Dcomp + MaskRow * (distPatche(:,:,1) - nDist(:,:,n)).^2;
    Dcomp(n) = norm(MaskRow *(distPatche(:,:,1) - nDist(:,:,n)));
  end

  %%% chosen
  for i=1:nNeighbors
  
    DcompBinaryList =  (Dcomp <= min(Dcomp)*(1+limitThreshold));
    [D, SimpleIdx] = min(Dcomp(DcompBinaryList));                   %comparação de quadrados (intensidade)
    SSDCoords = find(DcompBinaryList)                              %retorna o índice do valor diferente de 0 (O escolhido)     
    Chosendist(i) = SSDCoords(SimpleIdx);
    R(i) = D;
    Dcomp(Chosendist(i)) = 10000000;                                %Altera o valor minimo.
  end
  
  %--------------------------------------------------

end