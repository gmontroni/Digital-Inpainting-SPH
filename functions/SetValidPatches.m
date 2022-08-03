function [Patches, iCoordsVSample, jCoordsVSample] = ...
SetValidPatches(BinarySample, InputImage, winsize, halfwinsize) 

  % STEP 1.
  %---------------------------------------------------
  [iCoordsSample, jCoordsSample] = find(BinarySample);

  %%% Crop sample region to rectangle
  %%% Cortar a região da amostra para retângulo
  iRecMax = max(iCoordsSample)
  iRecMin = min(iCoordsSample)
  jRecMax = max(jCoordsSample)
  jRecMin = min(jCoordsSample)

  %%% Binary region containing binary external rectangle
  %%% Região binária contendo retângulo externo binário
  BinaryRectangle = BinarySample(iRecMin:iRecMax, jRecMin:jRecMax)
      
  %%% Get all binary window-patches into the big binary rectangle
  %%% Coloque todos os patches de janela binários no grande retângulo binário
  BinaryPatches = im2col(BinaryRectangle, [winsize winsize], 'sliding')
  [nPixelsByPatch, nPatches] = size(BinaryPatches)
  
  %%% Set list containing all valid region to a genuine candidate
  %for n = 1 : nPatches
  %%% Defina a lista contendo toda a região válida para um candidato genuíno
  % para n = 1: nPatches
  %  validList(n) = (sum(BinaryPatches(:,n)) == nPixelsByPatch);
  %end
  %validList = validList
  validList = (sum(BinaryPatches,1) == nPixelsByPatch)
  %---------------------------------------------------
  
  % STEP 1.1 (Update)
  %---------------------------------------------------
  %%% The use of heuristic when none of the candidates are found  
  %%% O uso da heurística quando nenhum dos candidatos é encontrado
  %%% This issue is addressed by taking the candidate with more 'ones'
  %%% Este problema foi solucionado levando o candidato com mais 'uns'
  if (sum(validList)==0)
     validList = (max(sum(BinaryPatches,1)) == sum(BinaryPatches,1))
  end
  %---------------------------------------------------
  
  % STEP 2.
  %---------------------------------------------------
  %%% Set number of rows and cols for which the patch can "walk on" the
  %%% big binary rectancle and get coords for pixels in the valid-sample
  %%% Defina o número de linhas e colunas pelas quais o patch pode "andar" no
  %%% grande retângulo binário e obtenha coordenadas para pixels na amostra válida
  walkPatches = size(BinaryRectangle) - winsize + 1
  [iCoordsVSample, jCoordsVSample] = ind2sub(walkPatches, find(validList))

  %%% Adjustment to correct the position w.r.t. big binary rectangle
  %%% Ajuste para corrigir a posição w.r.t. grande retângulo binário
  iCoordsVSample = iCoordsVSample + halfwinsize
  jCoordsVSample = jCoordsVSample + halfwinsize

  %%% Adjustment to correct the position w.r.t. all binary sample image
  %%% Ajuste para corrigir a posição w.r.t. toda imagem de amostra binária
  iCoordsVSample = iCoordsVSample + iRecMin - 1
  jCoordsVSample = jCoordsVSample + jRecMin - 1
  %---------------------------------------------------
  
  % STEP 3.
  %---------------------------------------------------
  %%% Get pixels into the input image with respect to big rectangle
  %%% Obtenha pixels na imagem de entrada em relação ao retângulo grande
  Rectangle = InputImage(iRecMin:iRecMax, jRecMin:jRecMax, :)
  ChansRec = size(Rectangle,3)

  %%% Get all window-patches sample into the rectangle and Get pixels from 
  %%% Obtenha todas as amostras de patches de janela no retângulo e obtenha pixels de
  %%% the valid window-patches with respect to input image
  %%% os patches de janela válidos em relação à imagem de entrada
  for n = 1 : ChansRec
    AuxPatches = im2col(Rectangle(:,:,n), [winsize winsize], 'sliding')
    Patches(:,:,n) = AuxPatches(:,validList)
  end
    %---------------------------------------------------

  % STEP 3.1 (Update)
  %---------------------------------------------
  %%% Heuristic when 'BinaryPatches' is empty (min(size(BinaryRectangle)) < winsize). 
  %%% Heurística quando 'BinaryPatches' estiver vazio (min (tamanho (BinaryRectangle)) <winsize).
  % This is solved by building an artificial patch with gaussian values 
  % Isso é resolvido através da criação de um patch artificial com valores gaussianos
  if isempty(BinaryPatches)
      clear Patches;
      for n = 1 : ChansRec
         ArtPatch = max(max(Rectangle(:,:,n)))*fspecial('gaussian', winsize, winsize);
         ArtPatch = im2col(ArtPatch, [winsize winsize], 'sliding');
         Patches(:,1,n) = ArtPatch;
         [sizeRows, sizeCols, ~] = size(InputImage);
         iCoordsVSample = round(sizeRows/2);
         jCoordsVSample = round(sizeCols/2);
      end    
  end
  %---------------------------------------------------
  
end 