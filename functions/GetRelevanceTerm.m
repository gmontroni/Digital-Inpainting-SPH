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