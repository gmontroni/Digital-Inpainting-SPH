function [SSDBinaryList, SSD, ChosenIdx] = ... 
SSDCompv2(Template, BinaryMask, gaussMask, Patches, nNeighbors)

  limitThreshold = 0.3;
  
  % STEP 1.
  %--------------------------------------------------
  %%% nPixelsPerPacth: Number of pixels of a candidate patch
  %%% nPatches:        Number of valid pixels from the sample area 
  [nPixelsPerPatch, nPatches, nChannels] = size(Patches);
  
  %%% The sum is computed only where 'BinaryMask' is true
  totalWeight = sum(sum(gaussMask(BinaryMask)));
  Mask = (gaussMask .* BinaryMask) / totalWeight;
  
  %%% Row-vector conversion
  MaskRow = Mask(:)';
  %--------------------------------------------------
  
  % STEP 2.
  %--------------------------------------------------
  for n = 1 : nChannels
    %%% 'Template' is replicated into 'TemplatePatches' [nPixelsPerPatch 
    %%% x nPatches x nChannels] to be compared with patches in 'Patches' 
    AuxTemplate(:,:,n) = reshape(Template(:,:,n), [nPixelsPerPatch 1]);
    TemplatePatches(:,:,n) = repmat(AuxTemplate(:,:,n), [1 nPatches]);  
  end
  %--------------------------------------------------

  % STEP 3.
  %--------------------------------------------------
  SSD = zeros(1, nPatches);
  %%% Computing SSD summing all image channels
  for n = 1 : nChannels
    SSD = SSD + MaskRow * (TemplatePatches(:,:,n) - Patches(:,:,n)).^2;   
  end
    
  %%% +1 to ensure that at least one of pixels from the sample will be
  %%% chosen
  for i=1:nNeighbors
  
    SSDBinaryList =  (SSD <= min(SSD)*(1+limitThreshold));
    [~, SimpleIdx] = min(SSD(SSDBinaryList));                   %comparação de quadrados (intensidade)
    SSDCoords = find(SSDBinaryList);                            %retorna o índice do valor diferente de 0 (O escolhido)     
    ChosenIdx(i) = SSDCoords(SimpleIdx);
    SSD(ChosenIdx(i)) = 10000000;                               %Altera o valor minimo.
  end
  
  %--------------------------------------------------

end