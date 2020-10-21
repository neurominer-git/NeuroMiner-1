function size_acc(y, x, init_size)

rand('seed', 0);
num_data = size(y,1);
if ~exist('init_size','var')
	init_size = round(num_data/4);
end
if num_data <= init_size
    disp('Data too small');
    exit;
end

if sum(y ~= floor(y)) > 0
    disp('Not a classification problem');
    exit;
end

if max(max(abs(x))) > 1
    disp('Warning: |x|_1 > 1; you may scale data for faster training');
end

perm_ind = randperm(num_data)';

xt_perm = x';
xt_perm = xt_perm(:,perm_ind);
y_perm = y(perm_ind);

% draw a new figure
figure; hold on; 
xlabel('Size of training subsets');
ylabel('Cross validation accuracy');
set(gca(), 'xminortick', 'on')

subset_size = init_size;
sizes = []; cv_accs = [];
while 1,
    x_subset = xt_perm(:, 1:subset_size)';
    y_subset = y_perm(1:subset_size);
    cv_accs = [cv_accs; train(y_subset, x_subset, '-v 5 -q')];
    sizes = [sizes; subset_size];
    semilogx(sizes, cv_accs, '*-');
    
    % set y-axis (accuracy) range at the first iteration
    if subset_size == init_size
        set(gca, 'ylim', [min(cv_accs(1)-5,90), 100]);
    end

    if subset_size >= num_data
        break;
    else
        subset_size = subset_size + 10;
        if subset_size>num_data,
            subset_size=num_data;
        end
    end
    drawnow;
end

hold off
