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

ITEM_RESULT = 0
ITEM_UP = 1
ITEM_NAME_UP = 2
ITEM_DOWN = 3
ITEM_NAME_DOWN = 4

def init_delta():
    return [0, 0, [], 0, []]

def modify_delta_up(deltas, locus, match_name):    
    if (not locus in deltas):
        deltas[locus] = init_delta()
    deltas[locus][ITEM_UP]  += 1;
    deltas[locus][ITEM_NAME_UP].append(match_name)

def modify_delta_down(deltas, locus, match_name):    
    if (not locus in deltas):
        deltas[locus] = init_delta()
    deltas[locus][ITEM_DOWN]  += 1;
    deltas[locus][ITEM_NAME_DOWN].append(match_name)

def calculate_deltas(row, required_chromosome, deltas):
    
    if (len(row) == NUMBER_OF_COLUMNS):
        match_name = row[COLUMN_MATCH_NAME]
        chromosome = int(row[COLUMN_CHROMOSOME])
        loc_beg = int(row[COLUMN_SEGMENT_BEGIN])
        loc_end = int(row[COLUMN_SEGMENT_END])
        if (chromosome != required_chromosome): # TODO: process all in one go
            return
        modify_delta_up(deltas, loc_beg, match_name)    
        modify_delta_down(deltas, loc_end, match_name)

def sum_deltas(deltas):
    number_of_matches_after_locus = 0
    for locus in sorted(deltas):
        number_of_matches_after_locus += deltas[locus][ITEM_UP] - deltas[locus][ITEM_DOWN]
        deltas[locus][ITEM_RESULT] = number_of_matches_after_locus

def plot_x_y(x, y, average):
    plt.xlabel("locus [Mb]")
    plt.ylabel("Number of matches")
    plt.title("Chromosome " + str(required_chromosome) + " average: " + str(int(average)))
    plt.plot(x, y)
    plt.show()


def plot_result(deltas, chromosome):
    x = []    
    y = []
    
    previous_locus = 0
    previous_result = 0
    
    locus_scaling_factor = 1000000
    total_no_of_samples = 0
    sum_of_samples = 0
    
    for locus in sorted(deltas):
        if (previous_locus > 0):    
            x.append(previous_locus / locus_scaling_factor)
            y.append(previous_result)
            x.append((locus - 1) / locus_scaling_factor)
            y.append(previous_result)
            no_of_samples = locus - previous_locus;
            total_no_of_samples += no_of_samples
            sum_of_samples += no_of_samples * previous_result
        previous_locus = locus
        previous_result = deltas[locus][ITEM_RESULT]
    
    plot_x_y(x, y, sum_of_samples / total_no_of_samples)

def print_name_list(name_list, header):
    if (len(name_list) > 0): 
        print (" ", header)
        for name in name_list:
            print("   ", name)


def print_result(deltas, chromosome, print_names):
    for locus in sorted(deltas):
        print (str(locus) + ";" + str(deltas[locus][ITEM_RESULT]))
        if (print_names):
            print_name_list(deltas[locus][ITEM_NAME_UP], "Begin of match:")
            print_name_list(deltas[locus][ITEM_NAME_DOWN], "End of match:")
            # print (deltas[locus][ITEM_NAME_UP], ";", deltas[locus][ITEM_NAME_DOWN])




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
parser.add_argument('--plot', dest='plot_result', action='store_true')
parser.add_argument('--print-names', dest='print_names', action='store_true')

args = parser.parse_args()

if (args.print_names and args.plot_result):
    print("Conflicting arguments --plot and --print-names", file=sys.stderr)
    exit() 

required_chromosome = int(args.required_chromosome)
deltas = dict()

if (args.csv_file_path):
    print("reading csv file", args.csv_file_path, file=sys.stderr)
    header = True
    with open(args.csv_file_path, newline='') as csvfile:
        csv_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in csv_reader:            
            if (header):
                header = False
                continue
            calculate_deltas(row, required_chromosome, deltas)

sum_deltas(deltas)

if (args.plot_result):
    plot_result(deltas, required_chromosome)
else:
    print_result(deltas, required_chromosome, args.print_names)


