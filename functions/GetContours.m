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