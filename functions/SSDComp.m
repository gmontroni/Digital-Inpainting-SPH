function [SSDBinaryList, SSD] = ... 
SSDComp(Template, BinaryMask, gaussMask, Patches)

  limitThreshold = 0.3;
  
  % STEP 1.
  %--------------------------------------------------
  %%% nPixelsPerPacth: Number of pixels of a candidate patch
  % nPixelsPerPacth: Número de pixels de um bloco candidato
  %%% nPatches:        Number of valid pixels from the sample area 
  % nPatches: Números de pixels válidos da area de amostra
  [nPixelsPerPatch, nPatches, nChannels] = size(Patches);
  
  %%% The sum is computed only where 'BinaryMask' is true
  % A soma é calculada apenas onde 'BinaryMask' é verdadeiro(BinaryMask
  % ==1)
  totalWeight = sum(sum(gaussMask(BinaryMask)));      %calcula o peso somente onde BinaryMask == 1
  Mask = (gaussMask .* BinaryMask) / totalWeight;     %aplica o filtro somente onde BinaryMask == 1 e divide pelo peso.
  
  %%% Row-vector conversion
  %Vetor Linha da Máscara.
  MaskRow = Mask(:)';
  %--------------------------------------------------
  
  % STEP 2.
  %--------------------------------------------------
  for n = 1 : nChannels
    %%% 'Template' is replicated into 'TemplatePatches' [nPixelsPerPatch 
    %%% x nPatches x nChannels] to be compared with patches in 'Patches' 
    AuxTemplate(:,:,n) = reshape(Template(:,:,n), [nPixelsPerPatch 1]);      %Template é colocado na forma de vetor coluna
    TemplatePatches(:,:,n) = repmat(AuxTemplate(:,:,n), [1 nPatches]);       %Concatena o Template para a comparação com os Patches(Todos os Blocos selecionados na região de amostra)
  end
  %--------------------------------------------------

  % STEP 3.
  %--------------------------------------------------
  SSD = zeros(1, nPatches);
  %%% Computing SSD summing all image channels
  for n = 1 : nChannels
    SSD = SSD + MaskRow * (TemplatePatches(:,:,n) - Patches(:,:,n)).^2;      %Calculo da métrica
  end
    
  %%% +1 to ensure that at least one of pixels from the sample will be
  %%% chosen
  SSDBinaryList =  (SSD <= min(SSD)*(1+limitThreshold)); %faz um teste lógico onde 0 é falso e 1 verdadeiro
  %--------------------------------------------------

end