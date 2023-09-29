#!/usr/bin/python3

import argparse
import csv
import matplotlib.pyplot as plt
import sys


def modify_delta(deltas, locus, direction):    
    if (locus in deltas):
        deltas[locus] += direction;
    else:
        deltas[locus] = direction;

def calculate_deltas(row, required_chromosome, deltas):
    
    direction_up = 1
    direction_down = -1

    if (len(row) == 10):
        chromosome = int(row[3])
        loc_beg = int(row[4])
        loc_end = int(row[5])
        if (chromosome != required_chromosome): # TODO: process all in one go
            return
        # print(loc_start, loc_end)
        modify_delta(deltas, loc_beg, direction_up)    
        modify_delta(deltas, loc_end, direction_down)

def sum_deltas(deltas, result):
    number_of_matches_after_locus = 0
    for locus in sorted(deltas):
        # print (locus, deltas[locus])
        number_of_matches_after_locus += deltas[locus]
        result[locus] = number_of_matches_after_locus

def plot_x_y(x, y, average):
    plt.xlabel("locus [MB]")
    plt.ylabel("Number of matches")
    plt.title("Chromosome " + str(required_chromosome) + " average: " + str(int(average)))
    plt.plot(x, y)
    plt.show()


def plot_result(result, chromosome):
    x = []    
    y = []
    
    previous_locus = 0
    previous_result = 0
    
    locus_scaling_factor = 1000000
    total_no_of_samples = 0
    sum_of_samples = 0
    
    for locus in sorted(result):
        if (previous_locus > 0):    
            x.append(previous_locus / locus_scaling_factor)
            y.append(previous_result)
            x.append((locus - 1) / locus_scaling_factor)
            y.append(previous_result)
            no_of_samples = locus - previous_locus;
            total_no_of_samples += no_of_samples
            sum_of_samples += no_of_samples * previous_result
        previous_locus = locus
        previous_result = result[locus]
    
    plot_x_y(x, y, sum_of_samples / total_no_of_samples)

def print_result(result, chromosome):
    for locus in sorted(result):
        print (str(locus) + ";" + str(result[locus]))




description='''
This program finds and shows pile-up regions in DNA test results. It analyses DNA match 
segments data provided as CSV file downloaded from MyHeritage and counts DNA matches
at particular loci. 
Result can be printed as table of locis and number of matches after the loci (default)
or plotted.
'''

parser = argparse.ArgumentParser(description=description)
parser.add_argument('--file', dest='csv_file_path')
parser.add_argument('--chromosome', dest='required_chromosome')
parser.add_argument('--plot', dest='plot_result', action='store_true')

args = parser.parse_args()

required_chromosome = int(args.required_chromosome)
result = dict()
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

sum_deltas(deltas, result)

if (args.plot_result):
    plot_result(result, required_chromosome)
else:
    print_result(result, required_chromosome)


