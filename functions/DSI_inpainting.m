
function [OutputImage, C, R] = DSI_inpainting(InputImage, Omega, win, winsize, winLarge, uband, numit)

%%% Description of the parameters
%-----------------------------------------------------------------------
% InputImage:  name of the image to be restored
% Omega:       binary inpainting region (1=to be filled, 0=other) 
% win:         window size to fill the target-exemplar 
% winsize:     window size to scan the sample (delimited by winLarge) 
% winLarge:    window size to determine the sampling region
% uband:       mapping image u (obtained from anisotropic diffusion)
% numit:       number of iterations used in the experiments
%-----------------------------------------------------------------------

%winsize: largura dos blocos.
%winlarge? largura da vizinhança.

warning off MATLAB:divideByZero


%%% Initialization
%-----------------------------------------------------------------------
%%% To properly work, this value must be integer
% Para funcionar corretamente, esse valor deve ser inteiro
halfwinsize = (winsize - 1) / 2;     %metade do win size.
nInpaint    = length(find(Omega));   %número de pixels a serem corrigidos.
nInpainted  = 0;                     %pixels corrigidos.
counter     = 0;                     %contador.

OutputImage = InputImage;            %imagem de entrada.
[sizeRows, sizeCols, ~] = size(InputImage);  %dimensão da imagem.
C = double(~Omega);
R = zeros([sizeRows sizeCols]);              %matriz de zeros com a dimensão da imagem de entrada.
indMat = reshape(1:winsize^2, winsize, winsize);  %matriz com elementos de 1 até winsize^2.

%%% Using gaussian distribuition to compute the weighted SSD metric
% Usando distribuição gaussiana para calcular a métrica SSD ponderada
gaussMask = fspecial('gaussian', winsize, winsize/6.4);      %suavizar a imagem.

%%% Computing the initial elements of the relevance term
% Computando os elementos iniciais do termo de relevância
[ulx, uly, ux, uy] = GetRelevanceTerm(uband);   %termo de relevância (pixel com prioridade)
%-----------------------------------------------------------------------



%==========================================================================
while (nInpainted < nInpaint) && (counter ~= numit)
  
  %%% Call function to identify the contour region to be filled
  %Chama a função para identificar a região de contorno de Omega
  %%% There are two choices: 0: Internal contour, 1: External contour
  %Existe duas chamadas; 0:Contorno interno e 1: Contorno externo
  %%% Notice that 'isConf' is a 'unitary' pixel contour (i.e.,length=1)
  isConf = GetContours(Omega,1);   %contorno
  
  %%% Computing the Confidence term C for all pixels from 'isConf'
  %Computando os termos de confiança de C para dos os pixels de isConf'
  %---------------------------------------------------------------
  for k = isConf'
    
    %%% Call function that gets a 'win x win' patch centered at 'k' 
    %Função de chamada que recebe um patch 'win x win' centralizado em 'k'
    Hk = GetPatch([sizeRows sizeCols], win, k);
      
    %%% Set indices from 'Hk' that properly aren't in the inpaint 
    % Definir índices de 'Hk' que não estão corretamente na pintura
    %%% domain. 'Hkfinal' is a vector
    Hkfinal = Hk(~(Omega(Hk)));
    
    %%% Compute the sum only in the pixels OUT of inpaint domain
    % Calcular a soma apenas nos pixels OUT do domínio inpaint
    %%% Notice that it's a recursive process for each C(k)
    C(k) = sum(C(Hkfinal)) / numel(Hk);    %equação 8 do paper ou eq 34 samara

  end 
  %---------------------------------------------------------------
  size(isConf)
  %%% Compute patch priorities: confidence * relevance 
  %Calcular prioridades do patch: confiança * relevância
  %%% ('Priorities' is a vector)
  R(isConf) = abs((ulx(isConf).* ux(isConf) + uly(isConf).* uy(isConf) + 0.001));   %equação 7 do paper ou eq 33 samara
  Priorities = C(isConf).*R(isConf);      %equação 9 do paper ou eq 32 samara (prioridades calculadas dos termos da borda)
    
  %%% Determine the local index with maximum priority
  % Determinar o índice local com prioridade máxima
  [~, ndx] = max(Priorities(:));                          %posição do pixel de maior prioridade da borda.
    
  %%% Compute the index of the pixel (p) that must be filled
  % Calcular o índice do pixel (p) que deve ser preenchido
  p = isConf(ndx(1));                                     %indice do termo com máxima prioridade.
  [iPixel, jPixel] = ind2sub([sizeRows sizeCols], p);     %localiza a posição (i,j) do pixel p na matrix InputImage (Imagem inicial).
  
  %%% Compute the i and j neighborhood w.r.t. the larger
  % Calcule o bairro i e j w.r.t. o maior
  %%% window, 'winsize', i.e. the window for pixel comparison
  % Janela, 'winsize', ou seja, a janela para comparação de pixels (janela
  % média, podemos deixar está janela para depois)
  iRange = max(iPixel-halfwinsize,1):min(iPixel+halfwinsize,sizeRows);       %tamanho da linha
  jRange = max(jPixel-halfwinsize,1):min(jPixel+halfwinsize,sizeCols);       %tamanho da coluna
  
  %%% Patch from the inpainted image and guidance mask (block of pixels)
  % Patch da imagem pintada e máscara de orientação (bloco de pixels)
  %podemos deixar está parte para depois, pois pertence ao winsize
  Template = OutputImage(iRange, jRange, :);    %Região escolhida para o retoque.
  BinaryMask = ~Omega(iRange, jRange);          %Máscara da região escolhida.
  
  %%% Determining the dynamic sampling region 'HpLarge' 
  % Determinando a região de amostragem dinâmica 'HpLarge'
  %%% (a 'winLarge x winLarge' patch centered at 'p') 
  % (um patch 'winLarge x winLarge' centralizado em 'p')
  HpLarge = GetPatch([sizeRows sizeCols], winLarge, p);  %Região de amostra (quadrado maior)
  OmegaSquare = Omega;
  OmegaSquare(HpLarge) = ones(size(HpLarge));            %Máscara do quadrado maior.
  DynamicSample = ~(OmegaSquare == Omega);               %onde as duas máscaras forem iguais, aplica o complementar (região da amostra na mascara).
  
  %%% Compute the candidate patches within the dynamic area 
  % Calcular os patches candidatos na área dinâmica
  %%% Note that Patches = [Rows: pixel values x Cols: candidate patches]
  % Observe que os patches = [Linhas: valores de pixel x Cols: patches candidatos]
  %%% 'iCoordsVSample' and 'jCoordsVSample': #cols of 'Patches' 
  % 'iCoordsVSample' e 'jCoordsVSample': #colunas de 'Patches'
  %%% and they store the central pixel of each col-patch of 'Patches'
  % e eles armazenam o pixel central de cada col-patch de 'Patches'
  
  [Patches, iCoordsVSample, jCoordsVSample] = ... 
  SetValidPatches(DynamicSample, OutputImage, winsize, halfwinsize);
  %função que calcula os patches. (imprimir)
  
  %figure, imshow(mat2gray(Patches));
  %size(Patches)
  %pause
  %%% Translating 'iRange' and 'jRange' to local system of coordinates
  % Traduzindo 'iRange' e 'jRange' para o sistema local de coordenadas
  iRangeValid = iRange-iPixel + halfwinsize+1;     
  jRangeValid = jRange-jPixel + halfwinsize+1;   
 
  %%% (a) Corner adjustmnet: gauss mask 'cut' w.r.t. fixed pacth 'Template'
  % (a) Rede de ajuste de canto: máscara gauss 'cut' w.r.t. pacth fixo 'Template'
  gaussM = gaussMask(iRangeValid, jRangeValid);    
    
  %%% (b) Corner adjustment: Patch 'cut' w.r.t. valid area of 'Template'
  % Ajuste de canto: Patch 'cut' w.r.t. área válida de 'Template'
  indTemplate = indMat(iRangeValid, jRangeValid);   %ajuste para a borda (indice do template)
  indTemplate = indTemplate(:);                     
  Patches = Patches(indTemplate, :,:);              %apenas os índices válidos.
  
  %%% Call function that calculates SSD distance between the 'fixed'
  % Função de chamada que calcula a distância SSD entre o 'fixo'
  %%% patch ('Template') and all the vec-col patches in 'Patches'
  [SSDBinaryList, SSD] = SSDComp(Template, BinaryMask, gaussM, Patches);            %calcula a distância dos blocos em relação ao bloco central
   
  %%% Determines the smallest pixel to be filled
  % Determina o menor pixel a ser preenchido
  [~, SimpleIdx] = min(SSD(SSDBinaryList));             %comparação de quadrados (intensidade)
  SSDCoords = find(SSDBinaryList);                      %retorna o índice do valor diferente de 0 (O escolhido)     
  ChosenIdx = SSDCoords(SimpleIdx);                     %Indice Escolhido
  
  [Hp, rows, cols] = GetPatch([sizeRows sizeCols], win, p);
  toFill = Omega(Hp);  % square matrix of 0 and 1
    
  %%% Search the best block on the sample
  % Procura o melhor bloco da amostra
  q = sub2ind([sizeRows sizeCols], iCoordsVSample(ChosenIdx), jCoordsVSample(ChosenIdx));   %Pega a submatriz e devolve o indice
  Hq = GetPatch([sizeRows sizeCols], win, q);
 
  %%% Propagate confidence and update. Here update C only
  % Propagar confiança e atualização. Aqui atualize apenas C
  %%% in the pixels that are containing into the region to be filled. 
  % nos pixels que estão contendo na região a ser preenchida.
  C(Hp(toFill)) = C(p);      %atualiza a matriz C
  
  %%% Update support relevance term     %atualiza os vetores (matriz R)
  % Atualizar termo de relevância do suporte% atualiza os vetores (matriz R)
  ulx(Hp(toFill)) = ulx(Hq(toFill));
  uly(Hp(toFill)) = uly(Hq(toFill));
  ux(Hp(toFill)) = ux(Hq(toFill));
  uy(Hp(toFill)) = uy(Hq(toFill)); 
    
  %%% Copy image data from Hq to Hp and
  % Copie os dados da imagem de Hq para Hp e
  %%% update fill region (region that is being filled)
  % atualizar região de preenchimento (região que está sendo preenchida)
  Omega(Hp(toFill)) = false;
  OutputImage(rows,cols,:) = UpdateImage(OutputImage,Hp,Hq,toFill);
  
  counter = counter + 1;
  
  nInpainted = nInpainted + length(find(toFill));
  fprintf('Pixels inpainted: %d / %d\n', nInpainted, nInpaint); 
  
  %figure, imshow(mat2gray(OutputImage))
  
  pause  
end %end wihle
%==========================================================================

end





%%% FUNCTIONS AND PROCEDURES (TO BE USED IN THE MAIN FUNCTION)
%==========================================================================
%==========================================================================

%%% Set input to the relevance term R
%=========================================================
function [ulx, uly, ux, uy] = GetRelevanceTerm(uband)
  
  %%% Computing lap(u)
  ulaux = zeros(size(uband));
  for k=1:size(uband,3)
     ulaux(:,:,k) = mat2gray(del2(uband(:,:,k)));
  end
  
  %%% Computing grad(lap(u))
  %ul = mean(ulaux,3);
  ul = max(ulaux,[],3);
  [ulx, uly] = gradient(ul);
  
  %%% Computing gradOrt(u) (with norm = 1)
  uax = zeros(size(uband));
  uay = zeros(size(uband));  
  for k=1:size(uband,3)
     [uax(:,:,k), uay(:,:,k)] = gradient(uband(:,:,k));
  end
  
  %ux = mean(mat2gray(uax),3);
  %uy = mean(mat2gray(uay),3);
  ux = max(mat2gray(uax),[],3);
  uy = max(mat2gray(uay),[],3);
   
  temp = ux;
  ux = -uy;
  uy = temp;

  norm = sqrt(ux.^2 + uy.^2);
  ux = ux./norm;
  uy = uy./norm;

end
%=========================================================
  
%%% Set valid patches from the sample
%=========================================================
function [Patches, iCoordsVSample, jCoordsVSample] = ...
SetValidPatches(BinarySample, InputImage, winsize, halfwinsize) 

  % STEP 1.
  %---------------------------------------------------
  [iCoordsSample, jCoordsSample] = find(BinarySample);

  %%% Crop sample region to rectangle
  iRecMax = max(iCoordsSample);
  iRecMin = min(iCoordsSample);
  jRecMax = max(jCoordsSample);
  jRecMin = min(jCoordsSample);

  %%% Binary region containing binary external rectangle
  BinaryRectangle = BinarySample(iRecMin:iRecMax, jRecMin:jRecMax);
      
  %%% Get all binary window-patches into the big binary rectangle
  BinaryPatches = im2col(BinaryRectangle, [winsize winsize], 'sliding');
  [nPixelsByPatch, nPatches] = size(BinaryPatches);
  
  %%% Set list containing all valid region to a genuine candidate
  %for n = 1 : nPatches
  %  validList(n) = (sum(BinaryPatches(:,n)) == nPixelsByPatch);
  %end
  %validList = validList
  validList = (sum(BinaryPatches,1) == nPixelsByPatch);
  %---------------------------------------------------
  
  % STEP 1.1 (Update)
  %---------------------------------------------------
  % The use of heuristic when none of the candidates are found  
  % This issue is addressed by taking the candidate with more 'ones'
  if (sum(validList)==0)
     validList = (max(sum(BinaryPatches,1)) == sum(BinaryPatches,1));
  end
  %---------------------------------------------------
  
  % STEP 2.
  %---------------------------------------------------
  %%% Set number of rows and cols for which the patch can "walk on" the 
  %%% big binary rectancle and get coords for pixels in the valid-sample
  walkPatches = size(BinaryRectangle) - winsize + 1;
  [iCoordsVSample, jCoordsVSample] = ind2sub(walkPatches, find(validList));

  %%% Adjustment to correct the position w.r.t. big binary rectangle
  iCoordsVSample = iCoordsVSample + halfwinsize;
  jCoordsVSample = jCoordsVSample + halfwinsize;

  %%% Adjustment to correct the position w.r.t. all binary sample image
  iCoordsVSample = iCoordsVSample + iRecMin - 1;
  jCoordsVSample = jCoordsVSample + jRecMin - 1;
  %---------------------------------------------------
  
  % STEP 3.
  %---------------------------------------------------
  %%% Get pixels into the input image with respect to big rectangle
  Rectangle = InputImage(iRecMin:iRecMax, jRecMin:jRecMax, :);
  ChansRec = size(Rectangle,3);

  %%% Get all window-patches sample into the rectangle and Get pixels from 
  %%% the valid window-patches with respect to input image
  for n = 1 : ChansRec
    AuxPatches = im2col(Rectangle(:,:,n), [winsize winsize], 'sliding');
    Patches(:,:,n) = AuxPatches(:,validList);
  end
  %---------------------------------------------------

  % STEP 3.1 (Update)
  %---------------------------------------------
  % Heuristic when 'BinaryPatches' is empty (min(size(BinaryRectangle)) < winsize). 
  % This is solved by building an artificial patch with gaussian values 
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
%=========================================================

%%% Get internal or external contour (not yet filled)
%=========================================================
function isConf = GetContours(Omega,IEOption)

  if (IEOption == 0)
  
    %%% Get internal border to be filled
    StrucElem = [0 1 0; 1 0 1; 0 1 0]; 
    BinaryContours = Omega - imerode(Omega,StrucElem);
  
  else
  
    %%% Get external border to be filled
    OmegaD = double(Omega);
    BinaryContours = conv2(OmegaD,[1,1,1;1,-8,1;1,1,1],'same')>0;
  
  end
  
  %%% Set simple indices of internal or external border  
  isConf = find(BinaryContours);    
  
end
%=========================================================
  
%%% Calculate SSD distance 
%=========================================================
function [SSDBinaryList, SSD] = ... 
SSDComp(Template, BinaryMask, gaussMask, Patches)

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
  SSDBinaryList =  (SSD <= min(SSD)*(1+limitThreshold));
  %--------------------------------------------------

end
%========================================================= 

% Returns the indices for a winf x winf patch centered at pixel k.
%========================================================= 
function [Hp, rows, cols] = GetPatch(sz, winf, pix)

  w = (winf - 1) / 2; 

  %%% Strategy to bypass the rest of the division above
  pix = pix-1;

  %%% floor: Rounds to the nearest integer 
  %%% Set the collumn position (y) of the k index 
  y = floor(pix/sz(1)) + 1; 
    
  %%% rem: Gets rest of the division (integer or real) 
  %%% Set the row position (x) of the k index 
  pix = rem(pix, sz(1));
  x = floor(pix) + 1;
  
  %%% Get range between (x,y)
  rows = max(x-w,1):min(x+w,sz(1));
  cols = (max(y-w,1):min(y+w,sz(2)))';

  %%% Set output block
  Hp = sub2ndx(rows, cols, sz(1));
  
end
%========================================================= 

% Converts the (rows,cols) subscript-style indices to Matlab index-style
% indices.  Unforunately, 'sub2ind' cannot be used for this.
%=========================================================
function Nt = sub2ndx(rows, cols, nTotalRows)

  %%% Repeat 'rows' length(cols) times. Same for the 'cols' 
  X = rows(ones(length(cols),1), :);
  Y = cols(:,ones(1,length(rows)));

  Nt = X+(Y-1)*nTotalRows;
  Nt = Nt';
  
end
%=========================================================

%%% Update OutputImage according to Hp(toFill) and Hq(toFill)
%=========================================================
function imUpdated = UpdateImage(LastVersionImage,Hp,Hq,toFill)

Hp(toFill) = Hq(toFill);
        
for n = 3:-1:1
  TempAux = LastVersionImage(:,:,n);
  imUpdated(:,:,n) = TempAux(Hp);
end

end
%=========================================================
