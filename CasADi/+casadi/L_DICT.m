function v = L_DICT()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = casadiMEX(0, 11);
  end
  v = vInitialized;
end
