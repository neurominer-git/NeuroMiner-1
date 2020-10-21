function [YTilde, YTeHat, YTeTilde] = myIsotonicRegression(YTr, YTe)

    startDir = pwd();

    ROOT_PATH = 'c:/users/aditya/documents/matlab/datasets/coil2000';
    LOOP_TYPE = 'roc';

    TRAIN_FILE = 'kdd08.train'; % 'coil-train.svml'
    TEST_FILE = 'kdd08.test'; % 'coil-test.svml'
    MODEL_FILE = 'kdd08.model'; % 'coil-model'
    
    TRAIN_FILE = sprintf('%s/%s', ROOT_PATH, TRAIN_FILE);
    TEST_FILE = sprintf('%s/%s', ROOT_PATH, TEST_FILE);
    MODEL_FILE = sprintf('%s/%s', ROOT_PATH, MODEL_FILE);

    %% Step 1: compute raw scores using ranking optimizer
    cd c:\cygwin\bin;    
    system(sprintf('%s/sofia-ml.exe --learner_type pegasos --loop_type %s --lambda 0.1 --training_file %s --model_out %s', ROOT_PATH, LOOP_TYPE, TRAIN_FILE, MODEL_FILE));
    system(sprintf('%s/sofia-ml.exe --model_in %s --test_file %s --results_file %s/results-train.txt', ROOT_PATH, MODEL_FILE, TRAIN_FILE, ROOT_PATH));
    system(sprintf('perl.exe %s/eval.pl %s/results-train.txt', ROOT_PATH, ROOT_PATH));
    
    cd(ROOT_PATH); %cd c:\users\aditya\documents\matlab\datasets\coil2000\;

    %% Step 2: sort training labels according to the scores
    x = load('results-train.txt');
    YHat = x(:,1);
    disp(sprintf('original train mse = %4.4f, auc = %4.4f, 0-1 = %4.4f', mean((YHat - (2*YTr-1)).^2), sampleError(YHat,YTr,'auc'), mean(((YHat - 0.5) .* YTr) >= 0)));

    [v,i] = sort(YHat);
    YTrSorted = YTr(i);

    %% Step 3: learn the isotonic mapping to this function
    YTilde = isotonic_regression(YTrSorted);
    
    % Add a slight perturbation to make sure there are no ties in the
    % isotonic function
    YTilde = YTilde + linspace(0,eps,length(YTilde))';
    
    disp(sprintf('isotonic train mse = %4.4f, auc = %4.4f', mean((YTilde - YTrSorted).^2), sampleError(YTilde,YTrSorted,'auc')));
    plot(1:length(YTilde), YTilde, 'o-', 1:length(YTrSorted), YTrSorted)
    
    %% Step 4: compute the raw scores on test examples
    cd c:\cygwin\bin;
    system(sprintf('%s/sofia-ml.exe --model_in %s --test_file %s --results_file %s/results-test.txt', ROOT_PATH, MODEL_FILE, TEST_FILE, ROOT_PATH));
    system(sprintf('perl.exe %s/eval.pl %s/results-test.txt', ROOT_PATH, ROOT_PATH));
    
    cd(ROOT_PATH); %cd c:\users\aditya\documents\matlab\datasets\coil2000\;    

    %% Step 5: replace each score with an isotonic counterpart
    x = load('results-test.txt');
    YTeHat = x(:,1);
    disp(sprintf('original test mse = %4.4f, auc = %4.4f, 0-1 = %4.4f', mean((YTeHat - (2*YTe - 1)).^2), sampleError(YTeHat,YTe,'auc'), mean(((YTeHat - 0.5) .* YTe) >= 0)));

    [v,i] = sort(YTeHat);
    YTeSorted = YTe(i);

    YTeTilde = zeros(size(YTe));
    for i=1:length(YTeHat)
        idx = bsearch(YTilde, v(i));
        YTeTilde(i) = YTilde(idx);
    end

    disp(sprintf('isotonic test mse = %4.4f, auc = %4.4f, 0-1 = %4.4f', mean((YTeTilde - YTeSorted).^2), sampleError(YTeTilde,YTeSorted,'auc'), mean(((YTeTilde - 0.5) .* YTeSorted) >= 0)));
    
    cd(startDir);
