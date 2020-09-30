import csv
import itertools
import os
import sys
import math
from math import log
import re
import random
import bisect


def read_seqs(path):
    container = []
    letters = {'A', 'C', 'G', 'T'}
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                append(line.strip().upper())
    return(container)


def read_peaks(path):
    container = []
    letters = {'A', 'C', 'G', 'T'}
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                append(line.strip().upper())
    return(container)


def read_seqs_with_complement(path):
    container = []
    letters = {'A', 'C', 'G', 'T'}
    append = container.append
    with open(path) as file:
        for line in file:
            if not line.startswith('>'):
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                append(line.strip().upper())
                append(complement(line.strip().upper()))
    return(container)


def complement(seq):
    return(seq.replace('A', 't').replace('T', 'a').replace('C', 'g').replace('G', 'c').upper()[::-1])


def read_fasta(path):
    fasta = list()
    letters = {'A', 'C', 'G', 'T'}
    with open(path, 'r') as file:
        for line in file:
            #print(line)
            if line.startswith('>'):
                line = line[1:].strip().split(':')
                record = dict()
                record['name'] = line[0]
                record['chromosome'] = line[2]
                coordinates_strand = line[3]

                start, end = re.findall(r'\d*-\d*', coordinates_strand)[0].split('-')
                record['start'] = start
                record['end'] = end

                strand = re.findall(r'\(.\)', coordinates_strand[:-3])
                if not strand == []:
                    record['strand'] = strand[0].strip('()')
                else:
                    record['strand'] = '+'
            else:
                line = ''.join([l if l in letters else 'N' for l in line.strip().upper()])
                record['seq'] = line.strip().upper()
                fasta.append(record)
    file.close()
    return(fasta)


def write_meme(output, tag, pfm, background, nsites):
    with open(output + '/' + tag + '.meme', 'w') as file:
        file.write('MEME version 4\n\nALPHABET= ACGT\n\nBackground letter frequencies\n')
        file.write('A {0} C {1} G {2} T {3}\n\n'.format(background['A'], background['C'],
                                                        background['G'], background['T']))
        file.write('MOTIF {0}\n'.format(tag))
        file.write(
            'letter-probability matrix: alength= 4 w= {0} nsites= {1}\n'.format(len(pfm['A']), nsites))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0:.8f}\t{1:.8f}\t{2:.8f}\t{3:.8f}\n'.format(i[0], i[1], i[2], i[3]))


def write_pwm(output, tag, pwm):
    with open(output + '/' + tag + '.pwm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pwm['A'], pwm['C'], pwm['G'], pwm['T']):
            file.write('{0}\t{1}\t{2}\t{3}\n'.format(i[0], i[1], i[2], i[3]))


def write_pfm(output, tag, pfm):
    with open(output + '/' + tag + '.pfm', 'w') as file:
        file.write('>{0}\n'.format(tag))
        for i in zip(pfm['A'], pfm['C'], pfm['G'], pfm['T']):
            file.write('{0:.9f}\t{1:.9f}\t{2:.9f}\t{3:.9f}\n'.format(i[0], i[1], i[2], i[3]))


def write_fasta(peaks, path):
    with open(path, 'w') as file:
        for index, p in enumerate(peaks):
            file.write('>{}\n'.format(index))
            file.write(p + '\n')
    return(0)


def creat_background(peaks, length_of_site, counter):
    shuffled_peaks = []
    number_of_sites = 0
    while counter > number_of_sites:
        peak = random.choice(peaks)
        shuffled_peak = ''.join(random.sample(peak, len(peak)))
        shuffled_peaks.append(shuffled_peak)
        number_of_sites += (len(''.join(shuffled_peak)) - length_of_site + 1) * 2
    return(shuffled_peaks)


def check_threshold_table(path):
    with open(path) as file:
        try:
            fpr = float(file.readline().strip().split()[1])
        except:
            fpr = -1
    return(fpr)
    

def write_scan(path, data):
    with open(path, 'w') as csvfile:
        fieldnames = ['chromosome', 'start', 'end', 'name', 'score', 'strand', 'site']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        #writer.writeheader()
        for line in data:
            writer.writerow(line)
    pass
    

# PWM MODEL
def make_pcm(motifs):
    matrix = {}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        matrix[i[0]] = []
    len_of_motif = len(motifs[0])
    for i in matrix.keys():
        matrix[i] = [0]*len_of_motif
    for i in range(len_of_motif):
        for l in motifs:
            matrix[l[i]][i] += 1
    return(matrix)


def make_pfm(pcm):
    number_of_sites = [0] * len(pcm['A'])
    for key in pcm.keys():
        for i in range(len(pcm[key])):
            number_of_sites[i] += pcm[key][i]
    pfm = dict()
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pfm[i[0]] = []
    first_key = list(pcm.keys())[0]
    nuc_pseudo = 1/len(pcm.keys())
    for i in range(len(pcm[first_key])):
        for nuc in pcm.keys():
            pfm[nuc].append((pcm[nuc][i] + nuc_pseudo) / (number_of_sites[i] + 1))
    return(pfm)


def make_pwm(pfm):
    pwm = {}
    background = {'A': 0.25,
                 'C': 0.25,
                 'G': 0.25,
                 'T': 0.25}
    mono_nucleotides = itertools.product('ACGT', repeat=1)
    for i in mono_nucleotides:
        pwm[i[0]] = []
    first_key = list(pfm.keys())[0]
    for i in range(len(pfm[first_key])):
        for j in pfm.keys():
            pwm[j].append(math.log(pfm[j][i] / background[j]))
    return(pwm)


def sites_to_pwm(sites):
    pcm = make_pcm(sites)
    pfm = make_pfm(pcm)
    pwm = make_pwm(pfm)
    return(pwm)


def read_pwm(path):
    with open(path, 'r') as file:
        inf = file.readline()
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
    file.close()
    return(pwm)