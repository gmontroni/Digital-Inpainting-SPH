function im = Kernel(LastVersionImage,winsize,h,R, nNeighbors,Hq,V)

im = zeros(winsize,winsize,3);
imAux = zeros(winsize,winsize,3);
[~,~,n] = size(LastVersionImage);

   for j=1:nNeighbors       
        W(j)= (1/(pi*(h^2)))*(exp(-((R(j))/(1))^2)) + 0.000001;
  end
  W
   Kernel = W*V;
   
  Kernel
    
%   for j=1:nNeighbors       
%       Kernel(j) = (1/(max(Kernel) - min(Kernel)))*(Kernel(j)- min(Kernel(:)))
%   end 
  %Kernel = mat2gray(Kernel);
  Kernel = normalize(Kernel);
  
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
  im
end