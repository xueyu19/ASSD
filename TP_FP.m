function [TP, FP] = TP_FP(h, h_pre)
    N = length(h);
    index_true = find(h);
    index_pre = find(h_pre);
    TP = length(intersect(index_true, index_pre));
    TPR = TP/length(index_true);
    FP = length(index_pre) - TP;
    FPR = FP/(N-length(index_true));
