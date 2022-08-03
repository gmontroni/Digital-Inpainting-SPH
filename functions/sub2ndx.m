function Nt = sub2ndx(rows, cols, nTotalRows)

  %%% Repeat 'rows' length(cols) times. Same for the 'cols' 
  X = rows(ones(length(cols),1), :);
  Y = cols(:,ones(1,length(rows)));

  Nt = X+(Y-1)*nTotalRows;
  Nt = Nt';
  
end