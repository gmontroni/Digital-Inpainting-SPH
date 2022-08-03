function im = DualKernel(LastVersionImage,winsize, nNeighbors,Hq,Ht)

alpha = 0.7;
imAux = zeros(winsize,winsize,3);
[~,~,n] = size(LastVersionImage);
  
  if n == 1
      imdist = zeros(winsize,winsize,n);
    for j=1:nNeighbors        
        TempAux = LastVersionImage(:,:,n);
        imAux(:,:,n) = TempAux(Hq(:,:,j));
        imdist(:,:,1) = imdist(:,:,1) + imAux(:,:,1)*Kernel(j);
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
    
       im = (alpha.*imdist + (1-alpha).*imint)./nNeighbors;
%      im = (alpha.*imdist + (1-alpha).*imint);
end