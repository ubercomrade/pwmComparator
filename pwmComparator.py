import os
import os.path
import sys
import subprocess
import argparse
import glob
import itertools
import shutil
from operator import itemgetter
from shutil import copyfile
from tools.get_threshold_for_pwm import get_threshold_for_pwm
from tools.scan_by_pwm import scan_by_pwm
from tools.prepare_bed import write_prepared_peaks
from tools.sites_intersection import sites_intersection
from tools.combine_results import combine_results
from tools.summary import write_peaks_classification
from tools.scan_best_by_pwm import scan_best_by_pwm
from tools.extract_sites import extract_sites
from tools.write_model import write_model
from lib.common import check_threshold_table, make_pwm

def prepare_data(path_to_genome, bed_path, bed, fasta):

    ########################
    #     GET TOP PEAKS    #
    ########################

    if not os.path.isfile(bed + '/' + 'peaks.bed'):
        #Get top training_sample_size bed peaks
        print('Prepare bed file')
        bed_out = bed + '/'
        write_prepared_peaks(bed_path, bed_out, 4, 'peaks')      
    else:
        print('{0} already exists'.format('peaks.bed'))

    ########################
    #     BED TO FASTA     #
    ########################

    if not os.path.isfile(fasta + '/' + 'peaks.fa'):
        #Bed peaks to fasta
        print('Get fasta from bed')
        bed_to_fasta(path_to_genome,
            bed + '/peaks.bed',
            fasta + '/peaks.fa')
    else:
        print('{0} already exists'.format('peaks.fa'))

    return(0)


def read_matrix(path):
    with open(path, 'r') as file:
        inf = file.readline()
        pwm = {'A': [], 'C': [], 'G': [], 'T': []}
        for line in file:
            line = line.strip().split('\t')
            for letter, value in zip(pwm.keys(), line):
                pwm[letter].append(float(value))
    file.close()
    return(pwm)


def calculate_thresholds_for_pwm(path_to_promoters, pwm_model_dir, thresholds_dir):
    if not os.path.isfile(thresholds_dir + '/pwm_model_thresholds.txt'):
        print('Calculate threshold for PWM based on promoters and fpr')
        get_threshold_for_pwm(path_to_promoters,
                pwm_model_dir + '/pwm_model.pwm',
                thresholds_dir + '/pwm_model_thresholds.txt')
    else:
        print('Thresholds for PWM already calculated')
    return(0)


def scan_peaks_by_pwm(fasta_test, model_path, model, scan, threshold_table_path, fpr):
    thr_pwm = get_threshold(threshold_table_path, fpr)
    pwm_scan_path = scan + '/{0}_{1:.2e}.bed'.format(model, fpr)
    print('Scan peaks by PWM with FPR: {0} THR: {1}'.format(fpr, thr_pwm))
    scan_by_pwm(fasta_test, model_path, thr_pwm, pwm_scan_path)
    return(0)


def bed_to_fasta(path_to_fa, path_to_bed, out):
    args = ['bedtools', 'getfasta' , '-s', '-name+',
            '-fi', path_to_fa,
            '-bed', path_to_bed,
            '-fo', out]
    r = subprocess.run(args, capture_output=True)
    pass


def get_threshold(path, fpr_for_thr):
    container = list()
    append = container.append
    with open(path, 'r') as file:
        for line in file:
            append(tuple(map(float, line.strip().split())))
    file.close()
    container = sorted(container, key=itemgetter(1))
    last_score, last_fpr = container[0]
    for line in container:
        if line[1] > fpr_for_thr:
            break
        else:
            last_score, last_fpr = line
    return(last_score)


def get_motif_length(models):
    with open(models + '/pwm_model/pwm_model.fasta', 'r') as file:
        for i in file:
            if i.startswith('>'):
                continue
            else:
                motif_length = len(i.strip())
                break
    file.close()
    return(motif_length)


def pipeline(models_names, models_paths, bed_path, fpr, \
    path_to_out, path_to_promoters, \
    path_to_genome, cpu_count):

    main_out = path_to_out
    cpu_count = cpu_count
    if not os.path.isdir(main_out):
        os.mkdir(main_out)
    models = main_out + '/models'
    thresholds = models + '/thresholds'
    fasta = main_out + '/fasta'
    bed = main_out + '/bed'
    scan = main_out + '/scan'
    scan_best = main_out + '/scan-best'
    results = main_out + '/results'
    
    ########################
    #      CREATE DIRS     #
    ########################

    if not os.path.isdir(models):
        os.mkdir(models)
    if not os.path.isdir(thresholds):
        os.mkdir(thresholds)
    if not os.path.isdir(fasta):
        os.mkdir(fasta)
    if not os.path.isdir(bed):
        os.mkdir(bed)
    if not os.path.isdir(scan):
        os.mkdir(scan)
    if not os.path.isdir(scan_best):
        os.mkdir(scan_best)
    if not os.path.isdir(results):
        os.mkdir(results)
    
    # PREPARE BED AND FASTA FILES #
    prepare_data(path_to_genome, bed_path, bed, fasta)

    peaks_fa = fasta + '/peaks.fa'
    peaks_bed = bed + '/peaks.bed'

    if models_names == 0 or len(models_names) != len(models_paths):
        models_names = ['model_{}'.format(i) for i in range(1, len(models_paths) + 1)]

    ### CALCULATE PWM MODEL ###
    for model, path in zip(models_names, models_paths):
        pwm_dir = models + '/{}_model/'.format(model)
        pwm_path = models + '/{0}_model/{0}_model.pwm'.format(model)
        pwm_threshold_table = thresholds + '/{}_model_thresholds.txt'.format(model)

        pfm = read_matrix(path)
        pwm = make_pwm(pfm)
        write_pwm(pwm_dir, model, pwm)
        
        # THRESHOLD
        if not os.path.isfile(thresholds_dir + '/pwm_model_thresholds.txt'):
            print('Calculate threshold for PWM based on promoters and fpr')
            get_threshold_for_pwm(path_to_promoters,
                    pwm_path,
                    pwm_threshold_table)
        else:
            print('Thresholds for {} already calculated'.format(model))
        
        check = check_threshold_table(pwm_threshold_table)
        if check < fpr:
            # SCAN
            scan_peaks_by_pwm(peaks_fa, pwm_model, model, scan, pwm_threshold_table, fpr)
            scan_best_by_pwm(scan_best + '/{}.scores.txt'.format(model),
                 pwm_model,
                 fasta_test)
        else:
            tools.remove(model)
            print('{} has poor table with thresholds'.format(model))
            print('GO to next model')
    ### END PWM ###

    # COMPARE SITES
    print('COMPARE SITES')
    pair_models = list(itertools.combinations(models, 2))
    for model1, model2 in pair_models:
        tag = 'compare'
        scan1 = scan + '/{0}_{1:.2e}.bed'.format(model1, fpr)
        scan2 = scan + '/{0}_{1:.2e}.bed'.format(model2, fpr)
        sites_intersection(bed_test, scan1, scan2, tag, model1, model2, results)

    # COMBINE SCAN
    list_bed_path = [scan + '/{0}_{1:.2e}.bed'.format(i, fpr) for i in tools]
    list_path_fpr_table = [thresholds + '/{}_model_thresholds.txt'.format(i) for i in models]
    combine_results(fasta_test, list_bed_path, list_path_fpr_table, tools, results + '/combined_scan.pro')

    # CALCULATE SUMMARY
    write_peaks_classification(results + '/combined_scan.pro', results + '/peaks_classification.tsv')
    print('Pipeline is finished!')
    tools = [t.upper() for t in tools]
    print('Results calculated for the next models:', *tools)
    

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('bed', action='store', help='path to BED file')
    parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10'], metavar='N',
         help='promoters of organism (hg38, mm10)')
    parser.add_argument('genome', action='store', help='path to genome fasta file')
    parser.add_argument('output', action='store', help='output dir')
    parser.add_argument('models', action='store', metavar='N', nargs='+',
         help='list of path to PFMs to use (./model1.pfm, ./model2.pfm ...)')
    parser.add_argument('-n' 'names', action='store', metavar='N', nargs='+', dest='names',
         help='list of names for models (name1, name2 ...)', required=False, default=[])
    parser.add_argument('-f', '--FPR', action='store', type=float, dest='fpr',
                        required=False, default=1.9*10**(-4), help='FPR, def=1.9*10^(-4)')
    parser.add_argument('-C', '--processes', action='store', type=int, dest='cpu_count',
                        required=False, default=4, help='Number of processes to use, default: 2')    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    bed_path = args.bed
    path_to_out = args.output
    fpr = args.fpr
    models_names = args.names
    models_paths = args.models

    organism = args.promoters
    path_to_genome = args.genome
    cpu_count = args.cpu_count

    this_dir, this_filename = os.path.split(__file__)
    if organism == 'mm10':
        path_to_promoters = os.path.join(this_dir, "promoters", "mm10.fasta")
    elif organism == 'hg38':
        path_to_promoters = os.path.join(this_dir, "promoters", "hg38.fasta")
    elif organism == 'tair10':
        path_to_promoters = os.path.join(this_dir, "promoters", "tair10.fasta")

    pipeline(models_names, models_paths, bed_path, fpr, path_to_out, path_to_promoters, \
        path_to_genome, cpu_count)

if __name__ == '__main__':
    main()
