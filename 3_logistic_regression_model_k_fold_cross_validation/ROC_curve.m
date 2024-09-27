function [TPR,FPR] = ROC_curve(variable,binary_class)
%ROC_curve: generates the true positive rate (TPR) and false postive rate
%(FPR) at varying thresholds values for a variable used to perform binary
%classification 

P = sum(binary_class); % all positives
N = length(binary_class) - P; % all negatives

steps = sort(variable); % sorts variable

for i = 1:length(steps)

    indx1 = find(variable<=steps(i)); % negatives
    indx2 = find(variable>steps(i)); % positives

    TP = sum(binary_class(indx2)); % true positive
    FP = length(indx2) - TP; % false positive
    
    TPR(i,:) = TP/P; % true positive rate
    FPR(i,:) = FP/N; % false positive rate

end

end