function do_binary_cross_validation(y, x, param, nr_fold)
len = length(y);
rand_ind = randperm(len);
dec_values = [];
labels = [];
for i = 1:nr_fold % Cross training : folding
  test_ind = rand_ind([floor((i-1)*len/nr_fold)+1:floor(i*len/nr_fold)]');
  train_ind = [1:len]';
  train_ind(test_ind) = [];
  model = svmtrain(y(train_ind),x(train_ind,:),param);
  [pred, acc, dec] = svmpredict(y(test_ind),x(test_ind,:),model);
  if model.Label(1) < 0;
    dec = dec * -1;
  end
  dec_values = vertcat(dec_values, dec);
  labels = vertcat(labels, y(test_ind));
end
disp(sprintf('Cross Validation:'));
validation_function(dec_values, labels);
