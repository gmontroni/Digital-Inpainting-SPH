function imUpdated = UpdateImage2(LastVersionImage,Hp,Hq,toFill)

 [m,l,p] = size(LastVersionImage);
 imUpdated = zeros(m,l,p);
        
for n = 3:-1:1
  TempAux = LastVersionImage(:,:,n)
  HqAux = Hq(:,:,n)
  TempAux(Hp(toFill)) = HqAux(toFill)
  %ver o que acontece nas linhas de cima e tentar arrumar embaixo
  imUpdated(:,:,n) = TempAux
end

end