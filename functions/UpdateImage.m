function imUpdated = UpdateImage(LastVersionImage,Hp,Hq,toFill)

Hp(toFill) = Hq(toFill);
        
for n = 3:-1:1
  TempAux = LastVersionImage(:,:,n);
  imUpdated(:,:,n) = TempAux(Hp);
end

end