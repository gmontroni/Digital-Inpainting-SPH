function [Hp, rows, cols] = GetPatch(sz, winf, pix)
  %Hk = GetPatch([sizeRows sizeCols], win, k)
  %k = isConf(1)
  w = (winf - 1) / 2; 

  %%% Strategy to bypass the rest of the division above
  %%% Estrat�gia para contornar o restante da divis�o acima
  pix = pix-1;

  %%% floor: Rounds to the nearest integer 
  %%% floor: Arredonda para o n�mero inteiro mais pr�ximo
  %%% Set the collumn position (y) of the k index
  %%% Define a posi��o da coluna (y) do �ndice k
  y = floor(pix/sz(1)) + 1    %localiza a coluna em que o indice k est�;
    
  %%% rem: Gets rest of the division (integer or real) 
  %%% rem: Obt�m o restante da divis�o (inteiro ou real)
  %%% Set the row position (x) of the k index 
  %%% Define a posi��o da linha (x) do �ndice k
  pix = rem(pix, sz(1))
  x = floor(pix) + 1          %localiza a linha em que o indice k est�;
  
  %%% Get range between (x,y)
  %%% Obter intervalo entre (x, y)
  
  rows = max(x-w,1):min(x+w,sz(1))        %prepara os elementos das linhas do quadrado Hp
  cols = (max(y-w,1):min(y+w,sz(2)))'     %prepara os elementos das colunas do quadrado Hp

  %%% Set output block
  %%% Definir bloco de sa�da
  Hp = sub2ndx(rows, cols, sz(1))         %Monta o quadrado Hp
  
end