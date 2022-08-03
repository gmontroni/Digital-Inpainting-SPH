function [] = teste(InputImage, Omega,win, winsize, winLarge, uband, numit)

%% Inicialização

%%% To properly work, this value must be integer
% Para funcionar corretamente, os valores precisam ser inteiros
halfwinsize = (winsize - 1) / 2;     %metade do win size.
nInpaint    = length(find(Omega));   %número de pixels a serem corrigidos.
nInpainted  = 0;                     %pixels corrigidos.
counter     = 0;                     %contador.

OutputImage = InputImage;            %imagem de entrada.
[sizeRows, sizeCols, ~] = size(InputImage);  %dimensão da imagem.
C = double(~Omega);                          %C é mascara invertida.
R = zeros([sizeRows sizeCols]);              %matriz de zeros com a dimensão da imagem de entrada.
indMat = reshape(1:winsize^2, winsize, winsize);  %matriz com elementos de 1 até winsize^2.

%%% Using gaussian distribuition to compute the weighted SSD metric
% Usando distribuição gaussiana para calcular a métrica SSD ponderada
gaussMask = fspecial('gaussian', winsize, winsize/6.4);      %suavizar a imagem.

%%% Computing the initial elements of the relevance term
% Computando os elementos iniciais do termo de relevância
%[ulx, uly, ux, uy] = GetRelevanceTerm(uband);   %termo de relevância (pixel com prioridade)

while (nInpainted < nInpaint) && (counter ~= numit)

  %%% Call function to identify the contour region to be filled
  %Chama a função para identificar a região de contorno de Omega
  %%% There are two choices: 0: Internal contour, 1: External contour
  %Existe duas chamadas; 0:Contorno interno e 1: Contorno externo
  %%% Notice that 'isConf' is a 'unitary' pixel contour (i.e.,length=1)
  isConf = GetContours(Omega,1);   %contorno
   
  for k = isConf'
    
    %%% Call function that gets a 'win x win' patch centered at 'k' 
    %Função de chamada que recebe um patch 'win x win' centralizado em 'k'
    Hk = GetPatch([sizeRows sizeCols], win, k);
      
    %%% Set indices from 'Hk' that properly aren't in the inpaint 
    % Definir índices de 'Hk' que não estão corretamente na pintura
    %%% domain. 'Hkfinal' is a vector
    Hkfinal = Hk(~(Omega(Hk)));  %quadrado sem o Omega.
    
    %%% Compute the sum only in the pixels OUT of inpaint domain
    % Calcular a soma apenas nos pixels OUT do domínio inpaint
    %%% Notice that it's a recursive process for each C(k)
    C(k) = sum(C(Hkfinal)) / numel(Hk);    %equação 8 do paper

  end   
    
counter = counter + 1;
  
nInpainted = nInpainted + length(find(toFill));
fprintf('Pixels inpainted: %d / %d\n', nInpainted, nInpaint);     
    
end

%% Funções

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

end