import random


def creat_table_bootstrap(true_scores, false_scores):
    cdef list table = []
    true_scores.sort(reverse=True)
    false_scores.sort(reverse=True)
    false_length = len(false_scores)
    true_length = len(true_scores)
    for tpr in [round(i * 0.01, 2) for i in range(5,105, 5)]:
        score = true_scores[round(true_length * tpr) - 1]
        actual_tpr = sum([1 if true_score >= score else 0 for true_score in true_scores]) / true_length
        fpr = sum([1 if false_score >= score else 0 for false_score in false_scores]) / false_length
        table.append({'Scores': score, 'TPR': tpr, 'ACTUAL_TPR': actual_tpr, 'FPR': fpr})
    return(table)


def score_pwm(str seq, dict pwm):
    cdef float score = 0
    cdef int position = 0 
    cdef str letter
    cdef int length = len(seq)
    for index in range(length):
        score += pwm[seq[index]][position]
        position += 1
    return(score)

    
def calculate_scores_pwm_bootstrap(list sites, dict pwm):
    cdef str site
    cdef list scores = []
    cdef int number_of_sites = len(sites)
    append = scores.append
    for i in range(number_of_sites):
        site = sites[i]
        scores.append(score_pwm(site, pwm))
    return(scores)


def calculate_scores_pwm_thresholds(list peaks, dict pwm, int length_of_site, float threshold):
    cdef str site
    cdef list scores = []
    cdef int i, N
    cdef int number_of_sites = 0
    cdef int number_of_peaks = len(peaks)
    append = scores.append
    for index in range(number_of_peaks):
        peak = peaks[index]
        N = len(peak) - length_of_site + 1
        for i in range(N):
            site = peak[i:length_of_site + i]
            if 'N' in site:
                continue
            number_of_sites += 1
            score = score_pwm(site, pwm)
            if score >= threshold:
                append(score)
    return(scores, number_of_sites)