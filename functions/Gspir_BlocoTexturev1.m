function [OutputImage,Omega,time,counter,w,alpha] = Gspir_BlocoTexturev1(InputImage, Omega, nNeighbors, winsize, winLarge, uband, numit)

warning off MATLAB:divideByZero

% Inicialização
% Para funcionar corretamente, esse valor deve ser inteiro
tic
halfwinsize = (winsize - 1) / 2;     %metade do win size.
nInpaint    = length(find(Omega));   %número de pixels a serem corrigidos.
nInpainted  = 0;                     %pixels corrigidos.
counter     = 0;                     %contador.

OutputImage = InputImage;                    %imagem de entrada.
[sizeRows, sizeCols, ~] = size(InputImage);  %dimensão da imagem.
C = double(~Omega);
D = zeros([sizeRows sizeCols]);              %matriz de zeros com a dimensão da imagem de entrada.
w = 0.1;
indMat = reshape(1:winsize^2, winsize, winsize);  %matriz com elementos de 1 até winsize^2.


% Usando distribuição gaussiana para calcular a métrica SSD ponderada
gaussMask = fspecial('gaussian', winsize, winsize/6.4);      %suavizar a imagem.


% Computando os elementos iniciais do termo de relevância
[ulx, uly, ux, uy] = GetRelevanceTerm(uband);   %termo de relevância (pixel com prioridade)

a = 1;
V = 1/nNeighbors;

figure(13), imshow(mat2gray(OutputImage));
pause
%-----------------------------------------------------------------------
while (nInpainted < nInpaint) && (counter ~= numit) %&& nInpainted <= 1000

    isConf = GetContours(Omega,1);   %contorno da região degradada.
    
%---------------------------------------------------------------
    for k = isConf'
    
    %Função de chamada que recebe um patch 'win x win' centralizado em 'k'
    Hk = GetPatch([sizeRows sizeCols], winsize, k);
      
    % Definir índices de 'Hk' que não estão corretamente na pintura
    Hkfinal = Hk(~(Omega(Hk)));
    
    % Calcular a soma apenas nos pixels OUT do domínio inpaint
    C(k) = sum(C(Hkfinal)) / numel(Hk);    %equação 8 do paper ou eq 34 samara

    end 
%---------------------------------------------------------------
 
  %Calcular prioridades do patch: confiança * relevância
  D(isConf) = abs((ulx(isConf).* ux(isConf) + uly(isConf).* uy(isConf) + 0.001));   %equação 7 do paper ou eq 33 samara
  R(isConf) = (1-w)*C(isConf) + w;
  Priorities = R(isConf)'.*D(isConf);      %equação 9 do paper ou eq 32 samara (prioridades calculadas dos termos da borda)
  %Priorities = isConf;
  
  % Determinar o índice local com prioridade máxima
  [~, ndx] = max(Priorities(:));                          %posição do pixel de maior prioridade da borda.
  
  % Calcular o índice do pixel (p) que deve ser preenchido
  p = isConf(ndx(1));                                     %indice do termo com máxima prioridade.
  [iPixel, jPixel] = ind2sub([sizeRows sizeCols], p);     %localiza a posição (i,j) do pixel p na matrix InputImage (Imagem inicial).
  
    %--------------------------Teste Visual----------------------------- 
     TESTE = OutputImage;  %teste visual
     TESTE6=TESTE; AUX=TESTE6(:,:,1); AUX(isConf)=0; TESTE6(:,:,1) = AUX; AUX=TESTE6(:,:,2); AUX(isConf)=0; TESTE6(:,:,2) = AUX; AUX=TESTE6(:,:,3); AUX(isConf)=0; TESTE6(:,:,3) = AUX;
     TESTE6(p) = 255;
     figure(7),imshow(uint8(TESTE6)); title('Priorities p'); %xlim([135 161]); ylim([206 232]);
     imwrite(uint8(TESTE6), 'Priorities.png')
  
  %--------------------------Teste Visual-----------------------------
  
  iRange = max(iPixel-halfwinsize,1):min(iPixel+halfwinsize,sizeRows);       %tamanho da linha
  jRange = max(jPixel-halfwinsize,1):min(jPixel+halfwinsize,sizeCols);       %tamanho da coluna
  
  iAuxDist = iRange'*ones(1,winsize);           %auxiliar para criar o vetor dist do Patche em i.
  jAuxDist = jRange'*ones(1,winsize);           %auxiliar para criar o vetor dist do Patche em j.
  distPatche(:,1) = iAuxDist(:); 
  distPatche(:,2) = sort(jAuxDist(:));          %vetor distância do Patche.
  
  Template = OutputImage(iRange, jRange, :);    %Região escolhida para o retoque.
  BinaryMask = ~Omega(iRange, jRange);          %Máscara da região escolhida.
  
  HpLarge = GetPatch([sizeRows sizeCols], winLarge, p);  %Região de amostra (quadrado maior)
  OmegaSquare = Omega;
  OmegaSquare(HpLarge) = ones(size(HpLarge));            
  DynamicSample = ~(OmegaSquare == Omega);               %Região da amostra na mascara.
  
  %função que calcula os patches. (imprimir)
  [Patches, iCoordsVSample, jCoordsVSample] = ... 
  SetValidPatches(DynamicSample, OutputImage, winsize, halfwinsize);

  iRangeValid = iRange-iPixel + halfwinsize+1;     
  jRangeValid = jRange-jPixel + halfwinsize+1;   
 
  gaussM = gaussMask(iRangeValid, jRangeValid);    
    
  indTemplate = indMat(iRangeValid, jRangeValid);   %indice do template.
  indTemplate = indTemplate(:);                     
  Patches = Patches(indTemplate, :,:);              %apenas os indices válidos.
  
  %calcula a distância dos blocos em relação ao bloco principal(intensidade) e seleciona os nNeighbors melhores.
  [SSDBinaryList, SSD, ChosenIndx] = SSDCompv2(Template, BinaryMask, gaussM, Patches, nNeighbors,winsize);            
  
  %calcula a distância dos blocos em relação ao bloco principal(coordenadas) e seleciona os nNeighbors melhores.
  [DcompBinaryList, r, Dcomp, ChosenIdx] = DistComp(distPatche, BinaryMask, Patches, sizeRows, sizeCols, winsize, iCoordsVSample, jCoordsVSample, halfwinsize, nNeighbors);

  h = a*max(r);
  BestPatchesDist = Patches(:,ChosenIdx,:);                            %seleciona os nNeighbors melhores patches pela distância.
  figure(10), imshow(mat2gray(BestPatchesDist)), title('Distância')
  BestPatchesInt = Patches(:,ChosenIndx,:);                             %seleciona os nNeighbors melhores patches pela intensidade.
  figure(11), imshow(mat2gray(BestPatchesInt)), title('Intensidade')
  figure(14), imshow(mat2gray(Patches(:,:,:))), title('Distância')
  figure(15), imshow(mat2gray(Patches(:,:,:))), title('Intensidade')
  %imwrite(mat2gray(BestPatchesInt), 'PatcheIntensidade.png')
  %imwrite(mat2gray(BestPatchesDist), 'PatcheDistancia.png')
  %imwrite(mat2gray(Patches(:,:,:)), 'Patches.png')
  
  %BestPatche = UpdatePatche(winsize,BestPatches,nNeighbors);          %média entre os nNeighbors patches.
  
  [Hp, rows, cols] = GetPatch([sizeRows sizeCols], winsize, p);    %monta o patche Hp.
  toFill = Omega(Hp);                                              %máscara onde os valores 'true' serão alterados.
  
  Hq = zeros(winsize,winsize,nNeighbors);
  for i=1:nNeighbors
  q(i) = sub2ind([sizeRows sizeCols], iCoordsVSample(ChosenIdx(i)), jCoordsVSample(ChosenIdx(i)));   %Pega a submatriz e devolve o indice
  Hq(:,:,i) = GetPatch([sizeRows sizeCols], winsize, q(i));
  end
  
  Ht = zeros(winsize,winsize,nNeighbors);
  for i=1:nNeighbors
  t(i) = sub2ind([sizeRows sizeCols], iCoordsVSample(ChosenIndx(i)), jCoordsVSample(ChosenIndx(i)));   %Pega a submatriz e devolve o indice
  Ht(:,:,i) = GetPatch([sizeRows sizeCols], winsize, t(i));
  end
%   h;
%   r;
%   nNeighbors;
%   V;

  %imwrite(OutputImage, 'problem.jpg')
%--------------------Matriz das Derivadas-----------------------------   
Hulx = DualKernel(ulx,winsize,nNeighbors,Hq,Ht);
Huly = DualKernel(uly,winsize,nNeighbors,Hq,Ht);
Hux = DualKernel(ux,winsize,nNeighbors,Hq,Ht);
Huy = DualKernel(uy,winsize,nNeighbors,Hq,Ht);

% Hulx = Kernel(ulx,winsize,h,r, nNeighbors,Hq,V);
% Huly = Kernel(uly,winsize,h,r, nNeighbors,Hq,V);
% Hux = Kernel(ux,winsize,h,r, nNeighbors,Hq,V);
% Huy = Kernel(uy,winsize,h,r, nNeighbors,Hq,V);

%---------------------------Kernel------------------------------------ 
[Hr,W] = Kernel(OutputImage,winsize,h,r, nNeighbors,Hq,V);
[Hr1,alpha] = DualKernel(OutputImage,winsize, nNeighbors,Hq,Ht);
%-------------------------------------------------------------------
  figure(20), imshow(uint8(Hr1))
  imwrite(uint8(Hr1(toFill)), 'BlocoResultante.png')
  % Propagar confiança e atualização. Aqui atualize apenas C
  % nos pixels que estão contendo na região a ser preenchida.
  C(Hp(toFill)) = C(p);      %atualiza a matriz C
  toFill
  
  % Atualizar termo de relevância do suporte% atualiza os vetores (matriz R)
  ulx(Hp(toFill)) = Hulx(toFill);
  uly(Hp(toFill)) = Huly(toFill);
  ux(Hp(toFill)) = Hux(toFill);
  uy(Hp(toFill)) = Huy(toFill);
%   alpha = 0.3*(ones(1,nNeighbors))
%   Hq = zeros(winsize^2,nNeighbors,3);
%   for i=1:nNeighbors
%   Hr1(:,:,i) = (alpha*(BestPatchesDist(:,:,i)') + (1-alpha)*BestPatchesInt(:,:,1)')/3;
%   end
 
  %figure(14), imshow((Hr1));
  
  % Copie os dados da imagem de Hq para Hp e
  % atualizar região de preenchimento (região que está sendo preenchida)
  Omega(Hp(toFill)) = false;
  OutputImage = UpdateImage2(OutputImage,Hp,Hr1,toFill);
  %OutputImage(rows,cols,:) = UpdateImage(I,Hp,Hq,toFill);
  
  counter = counter + 1;
  
  hold on
  nInpainted = nInpainted + length(find(toFill));
  fprintf('Pixels inpainted: %d / %d\n', nInpainted, nInpaint); 
  figure(13), imshow(mat2gray(OutputImage));
  %--------------------------Teste Visual----------------------------- 
   
     TESTE(Hp) = 155;
     TESTE1=TESTE; AUX=TESTE1(:,:,1); AUX(Hq(:,:,1))=255;  TESTE1(:,:,1)=AUX; AUX=TESTE1(:,:,2); AUX(Hq(:,:,1))=255; TESTE1(:,:,2)=AUX; AUX=TESTE1(:,:,3); AUX(Hq(:,:,1))=0; TESTE1(:,:,3)=AUX;
     TESTE2=TESTE; AUX=TESTE2(:,:,1); AUX(Hq(:,:,2))=255;  TESTE2(:,:,1)=AUX; AUX=TESTE2(:,:,2); AUX(Hq(:,:,2))=255; TESTE2(:,:,2)=AUX; AUX=TESTE2(:,:,3); AUX(Hq(:,:,2))=0; TESTE2(:,:,3)=AUX;
     TESTE3=TESTE; AUX=TESTE3(:,:,1); AUX(Hq(:,:,3))=255;  TESTE3(:,:,1)=AUX; AUX=TESTE3(:,:,2); AUX(Hq(:,:,3))=255; TESTE3(:,:,2)=AUX; AUX=TESTE3(:,:,3); AUX(Hq(:,:,3))=0; TESTE3(:,:,3)=AUX;
     TESTE4=TESTE; AUX=TESTE4(:,:,1); AUX(Hq(:,:,4))=255;  TESTE4(:,:,1)=AUX; AUX=TESTE4(:,:,2); AUX(Hq(:,:,4))=255; TESTE4(:,:,2)=AUX; AUX=TESTE4(:,:,3); AUX(Hq(:,:,4))=0; TESTE4(:,:,3)=AUX;
     TESTE5=TESTE; AUX=TESTE5(:,:,1); AUX(Hp)=0;  TESTE5(:,:,1)=AUX; AUX=TESTE5(:,:,2); AUX(Hp)=0; TESTE5(:,:,2)=AUX; AUX=TESTE5(:,:,3); AUX(Hp)=255; TESTE5(:,:,3)=AUX;
     TESTE7=TESTE; AUX=TESTE7(:,:,1); AUX(isConf)=0;       TESTE7(:,:,1)=AUX; AUX=TESTE7(:,:,2); AUX(isConf)=0;      TESTE7(:,:,2)=AUX; AUX=TESTE7(:,:,3); AUX(isConf)=0;    TESTE7(:,:,3)=AUX;
     
     figure(1), imshow(uint8(TESTE1)); title('nNeighbors Hq(1)'); %%xlim([135 161]); ylim([206 232])
     figure(2), imshow(uint8(TESTE2)); title('nNeighbors Hq(2)'); %xlim([135 161]); ylim([206 232]);
     figure(3), imshow(uint8(TESTE3)); title('nNeighbors Hq(3)'); %xlim([135 161]); ylim([206 232]);
     figure(4), imshow(uint8(TESTE4)); title('nNeighbors Hq(4)'); %xlim([135 161]); ylim([206 232]);
     figure(5), imshow(uint8(TESTE5)); title('Target Hp'); %xlim([135 161]); ylim([206 232]);
     imwrite(mat2gray(uint8(TESTE5)), 'Target.png')
    %--------------------------Teste Visual-----------------------------
    figure(6),imshow(uint8(OutputImage)); title('UpdateImage Hr'); %xlim([135 161]); ylim([206 232]);
    figure(8),imshow(uint8(TESTE7));      title('Contours'); %xlim([135 161]); ylim([206 232]);
    imwrite(mat2gray(OutputImage), 'UpdatedImage.png')
  
       pause

end %end wihle
time = toc;
imwrite(mat2gray(OutputImage), 'GSPIRBloco_output.png')
%==========================================================================
end %end function



%% Funções

%=========================================================
%% Set input to the relevance term R
% Configure a entrada para o termo de relevância R
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
  
%% Set valid patches from the sample
% Definir patches válidos da amostra
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

%% Get internal or external contour (not yet filled)
% Obter contorno interno ou externo (ainda não preenchido)
%=========================================================
function isConf = GetContours(Omega,IEOption)

  if (IEOption == 0)
  
    %%% Get internal border to be filled
    % Obter borda interna a ser preenchida
    StrucElem = [0 1 0; 1 0 1; 0 1 0]; 
    BinaryContours = Omega - imerode(Omega,StrucElem);
  
  else
  
    %%% Get external border to be filled
    % Obter borda externa a ser preenchida
    OmegaD = double(Omega);
    BinaryContours = conv2(OmegaD,[1,1,1;1,-8,1;1,1,1],'same')>0;
  
  end
  
  %%% Set simple indices of internal or external border  
  % Definir índices simples de borda interna ou externa
  isConf = find(BinaryContours);    
  
end
%=========================================================
  
%% Calculate SSD distance 
% Calcular a distância SSD
%=========================================================
function [SSDBinaryList, SSD, ChosenIdx] = ... 
SSDCompv2(Template, BinaryMask, gaussMask, Patches, nNeighbors,winsize)

  limitThreshold = 0.3;
  
  % STEP 1.
  %--------------------------------------------------
  %%% nPixelsPerPacth: Number of pixels of a candidate patch
  %%% nPatches:        Number of valid pixels from the sample area 
  [nPixelsPerPatch, nPatches, nChannels] = size(Patches);
  
  %%% The sum is computed only where 'BinaryMask' is true
  if winsize == 1
    totalWeight = 0.0558;
    %Mask = (gaussMask .* BinaryMask) / totalWeight;
    Mask = 1;
  else
    totalWeight = sum(sum(gaussMask(BinaryMask)));
    Mask = (gaussMask .* BinaryMask) / totalWeight;
  end
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
    [~, SimpleIdx] = min(SSD(SSDBinaryList));                %comparação de quadrados (intensidade)
    SSDCoords = find(SSDBinaryList);                         %retorna o índice do valor diferente de 0 (O escolhido)     
    ChosenIdx(i) = SSDCoords(SimpleIdx);
    SSD(ChosenIdx(i)) = 10000000;                               %Altera o valor minimo.
  end
  %--------------------------------------------------

end
%========================================================= 

%% Returns the indices for a winf x winf patch centered at pixel k.
% Retorna os índices para um patch winf x winf centralizado no pixel k.
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
  % Obter intervalo entre (x, y)
  rows = max(x-w,1):min(x+w,sz(1));
  cols = (max(y-w,1):min(y+w,sz(2)))';

  %%% Set output block
  Hp = sub2ndx(rows, cols, sz(1));
  
end
%========================================================= 

%% Converts the (rows,cols) subscript-style indices to Matlab index-style
% Converte os índices em estilo de subscrito (linhas, cols) em estilo de índice Matlab
% indices.  Unforunately, 'sub2ind' cannot be used for this.
% índices. Infelizmente, 'sub2ind' não pode ser usado para isso.
%=========================================================
function Nt = sub2ndx(rows, cols, nTotalRows)

  %%% Repeat 'rows' length(cols) times. Same for the 'cols' 
  X = rows(ones(length(cols),1), :);
  Y = cols(:,ones(1,length(rows)));

  Nt = X+(Y-1)*nTotalRows;
  Nt = Nt';
  
end
%=========================================================

%% Update OutputImage according to Hp(toFill) and Hq(toFill)
% Atualize OutputImage de acordo com Hp (toFill) e Hq (toFill)
%=========================================================
function imUpdated = UpdateImage(LastVersionImage,Hp,Hq,toFill)

Hp(toFill) = Hq(toFill);
        
for n = 3:-1:1
  TempAux = LastVersionImage(:,:,n);
  imUpdated(:,:,n) = TempAux(Hp);
end

end

%% Update OutputImage according to Hp(toFill) and Hq(toFill)
% Atualize OutputImage de acordo com Hp (toFill) e Hq (toFill) for distComp
%=========================================================
function imUpdated = UpdateImage2(LastVersionImage,Hp,Hq,toFill)

 [m,l,p] = size(LastVersionImage);
 imUpdated = zeros(m,l,p);
        
for n = 3:-1:1
  TempAux = LastVersionImage(:,:,n);
  HqAux = Hq(:,:,n);
  TempAux(Hp(toFill)) = HqAux(toFill);
  imUpdated(:,:,n) = TempAux;
end

end



%=========================================================

%=========================================================

function [DcompBinaryList, R, Dcomp, Chosendist] = ... 
DistComp(distPatche, BinaryMask, Patches, sizeRows, sizeCols, winsize, iCoordsVSample, jCoordsVSample, halfwinsize, nNeighbors)

  [nPixelsPerPatch, nPatches, nChannels] = size(Patches);
  limitThreshold = 0.3;
  MaskRow = BinaryMask(:)';
  nDist = zeros(nPixelsPerPatch,2,nPatches);
  
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
  %%% Computing SSD summing all image channels
  for n = 1 : nPatches
    %Dcomp = Dcomp + MaskRow * (distPatche(:,:,1) - nDist(:,:,n)).^2;
    Dcomp(n) = norm(MaskRow *(distPatche(:,:,1) - nDist(:,:,n)));
  end
    
  %%% chosen
  for i=1:nNeighbors
  
    DcompBinaryList =  (Dcomp <= min(Dcomp)*(1+limitThreshold));
    [D, SimpleIdx] = min(Dcomp(DcompBinaryList));                   %comparação de quadrados (intensidade)
    SSDCoords = find(DcompBinaryList);                              %retorna o índice do valor diferente de 0 (O escolhido)     
    Chosendist(i) = SSDCoords(SimpleIdx);
    R(i) = D;
    Dcomp(Chosendist(i)) = 10000000;                                %Altera o valor minimo.
  end
  
  %--------------------------------------------------

end
%--------------------------------------------------
function BestPatche = UpdatePatche(winsize,BestPatches,nNeighbors)

BestPatche = zeros(winsize^2,1,3);
for n = 1:nNeighbors
    
  BestPatche = BestPatche + BestPatches(:,n,:);
end

BestPatche = BestPatche/nNeighbors;

end
%--------------------------------------------------

function [im,W] = Kernel(LastVersionImage,winsize,h,R, nNeighbors,Hq,V)

im = zeros(winsize,winsize,3);
imAux = zeros(winsize,winsize,3);
[~,~,n] = size(LastVersionImage);

   for j=1:nNeighbors       
        W(j)= (1/(pi*(h^2)))*(exp(-((R(j))/(1))^2));
   end
   
%    W = max(0.5, W);
   Kernel = W*V;
    
%   for j=1:nNeighbors       
%       Kernel(j) = (1/(max(Kernel) - min(Kernel)))*(Kernel(j)- min(Kernel(:)))
%   end 
  %Kernel = mat2gray(Kernel);
   Kernel = normalize(Kernel,'norm',1);
   W = Kernel;
  if n == 1
      im = zeros(winsize,winsize,n);
    for j=1:nNeighbors        
        TempAux = LastVersionImage(:,:,n);
        imAux(:,:,n) = TempAux(Hq(:,:,j));
        im(:,:,1) = im(:,:,1) + imAux(:,:,1)*Kernel(j);
    end
  else
    for j=1:nNeighbors        
        for n = 3:-1:1
            TempAux = LastVersionImage(:,:,n);
            imAux(:,:,n) = TempAux(Hq(:,:,j));
        end
        im(:,:,:) = im(:,:,:) + imAux(:,:,:)*Kernel(j);
    end
  end
%   im
end

function [im,alpha] = DualKernel(LastVersionImage,winsize, nNeighbors,Hq,Ht)

alpha = 0.3;
imAux = zeros(winsize,winsize,3);
[~,~,n] = size(LastVersionImage);
  
  if n == 1
      imdist = zeros(winsize,winsize,n);
    for j=1:nNeighbors        
        TempAux = LastVersionImage(:,:,n);
        imAux(:,:,n) = TempAux(Hq(:,:,j));
        imdist(:,:,1) = imdist(:,:,1) + imAux(:,:,1);
    end
    imint = zeros(winsize,winsize,n);
     for j=1:nNeighbors        
         TempAux = LastVersionImage(:,:,n);
         imAux(:,:,n) = TempAux(Ht(:,:,j));
         imint(:,:,1) = imint(:,:,1) + imAux(:,:,1);
     end
  else
      imdist = zeros(winsize,winsize,3);
      imint = zeros(winsize,winsize,3);
    for j=1:nNeighbors        
        for n = 3:-1:1
            TempAux = LastVersionImage(:,:,n);
            imAux(:,:,n) = TempAux(Hq(:,:,j));
        end
        imdist(:,:,:) = imdist(:,:,:) + imAux(:,:,:);
    end
        
    for j=1:nNeighbors        
        for n = 3:-1:1
            TempAux = LastVersionImage(:,:,n);
            imAux(:,:,n) = TempAux(Ht(:,:,j));
        end
        imint(:,:,:) = imint(:,:,:) + imAux(:,:,:);
    end
   
  end
     
      im = (alpha*imdist + (1-alpha)*imint)/nNeighbors;
      %im = (alpha.*imdist2 + (1-alpha).*imint)/nNeighbors;
%     im = (alpha.*imdist + (1-alpha).*imint);
end