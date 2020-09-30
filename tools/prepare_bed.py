import os
import random


def read_file(path):
    peaks = []
    with open(path, 'r') as file:
        for line in file:
            if line.isspace():
                continue
            else:
                line = line.strip().split()
                peaks.append(line)
    file.close()
    return(peaks)


def clear_peaks(peaks):
    peaks = [i for i in peaks if len(i[0]) < 6 and i[0] != "chrMT"]
    return(peaks)


def prepare_peaks_bed6(peaks, col):
    for index, line in enumerate(peaks):
        line[3] = 'peaks_' + str(index)
    return(peaks)


def prepare_peaks_bed3(peaks):
    key = len(peaks[0]) < 4
    for index, line in enumerate(peaks):
        if key:
            line.append('peaks_' + str(index))
        else:
            line[3] = 'peaks_' + str(index)
    return(peaks)


def get_legths(data):
    l = list()
    for line in data:
        l.append(int(line[2]) - int(line[1]))
    return(l)



def write_length(path, data):
    with open(path, "w") as file:
        for line in data:
            file.write("{0}\n".format(line))
    file.close()
    return(0)


def write_prepared_peaks(path, output, col, tag):
    if not os.path.isdir(output):
        os.mkdir(output)
    peaks = read_file(path)
    peaks = clear_peaks(peaks)
    if len(peaks[0]) >= 5:
        peaks = prepare_peaks_bed6(peaks, col)
    else:
        peaks = prepare_peaks_bed3(peaks)
    l = get_legths(peaks)
    write_length(output + '/' + tag + '.length.txt', l)
    with open(output + '/' + tag + '.bed', 'w') as file:
        for i in peaks:
            i[-1] = str(i[-1])
            file.write('\t'.join(i) + '\n')
    file.close()
    return(0)