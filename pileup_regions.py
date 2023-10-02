#!/usr/bin/python3

import argparse
import csv
import matplotlib.pyplot as plt
import sys

NUMBER_OF_COLUMNS = 10 # Number of columns in input CSV
COLUMN_CHROMOSOME = 3 
COLUMN_SEGMENT_BEGIN = 4
COLUMN_SEGMENT_END = 5
COLUMN_MATCH_NAME = 2
COLUMN_CM = 8

ITEM_NUMBER_OF_MATCHES = 0
ITEM_SEGMENT_BEGIN = 1
ITEM_SEGMENT_BEGIN_NAMES = 2
ITEM_SEGMENT_END = 3
ITEM_SEGMENT_END_NAMES = 4

def init_location(locations, bp):
    if (not bp in locations):
        locations[bp] = [0, 0, [], 0, []]

def register_segment_begin(locations, bp, match_name):    
    init_location(locations, bp)
    locations[bp][ITEM_SEGMENT_BEGIN] += 1
    locations[bp][ITEM_SEGMENT_BEGIN_NAMES].append(match_name)

def register_segment_end(locations, bp, match_name):    
    init_location(locations, bp)
    locations[bp][ITEM_SEGMENT_END] += 1
    locations[bp][ITEM_SEGMENT_END_NAMES].append(match_name)

def register_segment(row, required_chromosome, cm_limit, locations):
    if (len(row) != NUMBER_OF_COLUMNS):
        return
    match_name = row[COLUMN_MATCH_NAME]
    chromosome = int(row[COLUMN_CHROMOSOME])
    centimorgans = float(row[COLUMN_CM])
    loc_beg = int(row[COLUMN_SEGMENT_BEGIN])
    loc_end = int(row[COLUMN_SEGMENT_END])
    if (chromosome != required_chromosome): # TODO: process all in one go
        return
    if (centimorgans < cm_limit):
        return
    register_segment_begin(locations, loc_beg, match_name)    
    register_segment_end(locations, loc_end, match_name)

def calculate_number_of_matches(locations):
    number_of_matches_after_locus = 0
    for bp in sorted(locations):
        location = locations[bp]
        number_of_matches_after_locus += location[ITEM_SEGMENT_BEGIN] - location[ITEM_SEGMENT_END]
        location[ITEM_NUMBER_OF_MATCHES] = number_of_matches_after_locus

def plot_x_y(x, y, average):
    plt.xlabel("locus [Mb]")
    plt.ylabel("Number of matches")
    plt.title("Chromosome " + str(required_chromosome) + " average: " + str(int(average)))
    plt.plot(x, y)


def plot_matches(locations, chromosome):
    x = []    
    y = []
    
    previous_bp = 0
    previous_matches = 0
    
    bp_scaling_factor = 1000000
    total_no_of_samples = 0
    sum_of_samples = 0
    
    for bp in sorted(locations):
        if (previous_bp > 0):    
            x.append(previous_bp / bp_scaling_factor)
            y.append(previous_matches)
            x.append((bp - 1) / bp_scaling_factor)
            y.append(previous_matches)
            no_of_samples = bp - previous_bp;
            total_no_of_samples += no_of_samples
            sum_of_samples += no_of_samples * previous_matches
        previous_bp = bp
        previous_matches = locations[bp][ITEM_NUMBER_OF_MATCHES]
    
    plot_x_y(x, y, sum_of_samples / total_no_of_samples)

def print_name_list(name_list, header):
    if (len(name_list) > 0): 
        print (" ", header)
        for name in name_list:
            print("   ", name)


def print_matches(locations, chromosome, print_names):
    for bp in sorted(locations):
        print (str(chromosome) + ";" + str(bp) + ";" + str(locations[bp][ITEM_NUMBER_OF_MATCHES]))
        if (print_names):
            print_name_list(locations[bp][ITEM_SEGMENT_BEGIN_NAMES], "Begin of match:")
            print_name_list(locations[bp][ITEM_SEGMENT_END_NAMES], "End of match:")
            # print (locations[bp][ITEM_SEGMENT_BEGIN_NAMES], ";", locations[bp][ITEM_SEGMENT_END_NAMES])




description='''
This program finds and shows pile-up regions in DNA test results. It analyses DNA match 
segments data provided as CSV file downloaded from MyHeritage and counts DNA matches
at particular loci. 
Result can be printed as table of locis and number of matches after the locus (default)
or plotted.
'''


parser = argparse.ArgumentParser(description=description)
parser.add_argument('--file', dest='csv_file_path')
parser.add_argument('--chromosome', dest='required_chromosome')
parser.add_argument('--cm-limit', dest='cm_limit')
parser.add_argument('--plot', dest='plot_matches', action='store_true')
parser.add_argument('--print-names', dest='print_names', action='store_true')
parser.add_argument('--save-picture', dest='save_picture')

args = parser.parse_args()

if (args.print_names and args.plot_matches):
    print("Conflicting arguments --plot and --print-names", file=sys.stderr)
    exit() 

cm_limit = 0
if (args.cm_limit):
    cm_limit = float(args.cm_limit)


required_chromosome = int(args.required_chromosome)
locations = dict()

if (args.csv_file_path):
    print("reading csv file", args.csv_file_path, file=sys.stderr)
    header = True
    with open(args.csv_file_path, newline='', encoding='utf-8') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in csv_reader:            
            if (header):
                header = False
                continue
            register_segment(row, required_chromosome, cm_limit, locations)

calculate_number_of_matches(locations)

if (args.plot_matches or args.save_picture):
    plot_matches(locations, required_chromosome)
    if (args.save_picture):
        plt.savefig(args.save_picture)
    else:
        plt.show()
else:
    print_matches(locations, required_chromosome, args.print_names)


