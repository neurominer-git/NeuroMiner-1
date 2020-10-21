function [pred ret dec] = do_binary_predict(y, x, model)
[pred acc dec] = svmpredict(y, x, model);
if model.Label(1) < 0;
  dec = dec * -1;
end
ret = validation_function(dec, y);
