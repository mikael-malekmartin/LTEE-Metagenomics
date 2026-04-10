#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
import textwrap
from contextlib import nullcontext

import utils_LTEE

cwd = os.getcwd() # Get current working directory to save plot and datas in the output folder


# Create an argument parser to get the population name, gene name and options for the plot and datas saving
parser = argparse.ArgumentParser(description=textwrap.dedent('''\
                                                Generate a plot of mutation for a specific population and a specific gene.
                                                Default output folder for datas and plots : [APP DIR]/output/
                                                '''),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-p", "--population", help="Population name. Format : m/p1-6 (i.e. m1 or p6)", required=True)
requiredNamed.add_argument("-g", "--gene", help="The name of a gene (i.e. topA)", required=True)
parser.add_argument("-i", "--include_intergenic", help="[Optionnal] Include intergenic mutation before or after gene of interest", action="store_true")
parser.add_argument("-s", "--show_plot", action="store_true", help="[Optionnal] Show the plot, to edit before save it (specified destination)")
parser.add_argument("-d", "--datas_save", action="store_true", help="[Optionnal] Save datas as csv coma separate in default output folder")
parser.add_argument("--xkcd", action="store_true", help="[Optionnal] Comic plot")

args = parser.parse_args()

config = vars(args)

pop = config.get('population').casefold()
gene = config.get('gene')
plot_show = config.get('show_plot')
include_intergenic = config.get('include_intergenic')
datas_save = config.get('datas_save')
xkcd = config.get('xkcd')

os.makedirs('output', exist_ok=True) # Create folders if not exist for output (plot and datas)

txtFile = f'input/{pop}_annotated_timecourse.txt'


# Check if the given population name is correct (m/p1-6) and raise an error if not correct
pop_list = ['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 
            'p1', 'p2', 'p3', 'p4', 'p5', 'p6']

if pop not in pop_list:
    raise ValueError(f"'{pop}' is not a correct population name. Try with the following format m/p1-6 (i.e. m1 or p6)")

# Check if the given file exist in the input folder and raise an error if not exist
try :
    open(txtFile, 'r')
except FileNotFoundError as exc:
    raise FileNotFoundError(
        f"Error: No file '{txtFile}' found. Check if the file is in the input folder."
    ) from exc

# Convert the given argument 'pop' to Population name
minusOrPlus = pop[0]

if minusOrPlus == 'm' :
    minusOrPlus = '-'
elif minusOrPlus == 'p' :
    minusOrPlus = '+'

popName = f'Ara{minusOrPlus}{pop[-1]}'

# Import datas from file to pandas dataframe
df = pd.read_csv(txtFile, delimiter=',')
df_columns = df.columns.to_list()

# Add columns to dataframe
if include_intergenic == True:
    insertCol = ['Intergenic between', 'AA change']
else: 
    insertCol = ['AA change']

df = pd.concat([df, pd.DataFrame(columns=insertCol)], axis=1)

df = df[1:] # Remove index 0 because useless

# Add information about intergenic sequences (Gene before and after)
if include_intergenic == True:
    intergenic_range = []
    intergenic = False
    geneBefore = ''
    geneAfter = ''

    for index, row in df.iterrows():
        if row[' Gene'] == ' intergenic' and index != 0:
            intergenic_range.append(index)
            if intergenic != True :
                geneBefore = df.loc[index - 1, ' Gene']
                intergenic = True
        else:
            geneAfter = df.loc[index, ' Gene']
            intergenic = False
            for i in intergenic_range:
                df.at[i, 'Intergenic between'] = f'{geneBefore}-{geneAfter}'
            intergenic_range, geneBefore, geneAfer = [], '', ''

if include_intergenic == True:
    df_filtered = df[(df[' Gene'].str.contains(gene, case=False)) | (df['Intergenic between'].str.contains(gene, case=False))]
else :
    df_filtered = df[df[' Gene'].str.contains(gene, case=False)]

# Raise an error if no mutation found or incorrect gene name 
if len(df_filtered) == 0:
    raise ValueError(f'Error: No mutation in {gene} found in population {popName} or incorrect gene name')

# Define range
range_timecourse = [i for i in range(0,60001,500)]

# Create a new dataframe for plot with range as index and mutation frequency as columns
df_plot = pd.DataFrame(index = range_timecourse)

# Get information about gene of interest from genbank file and add AA 
genbank606 = utils_LTEE.genbank('input/REL606.6.gbk')
cds606 = utils_LTEE.cds_features(genbank606)
geneList, seq_ref, ref_codons, strand = utils_LTEE.geneProfile(genbank606, cds606, gene)

# Calculate mutation frequency for each mutation and add it to the plot dataframe with a header containing position, allele, annotation and AA change if missense mutation
for index, row in df_filtered.iterrows():
    frequency_list = []

    for i in range_timecourse:
        try :
            ac = row[f' AC:{i}']
            dp = row[f' DP:{i}']
        except KeyError:
            row[f' AC:{i}'] = np.nan
            row[f' DP:{i}'] = np.nan
        if dp != 0:
            frequency = ac / dp
        else:
            frequency = np.nan
 
        frequency_list.append(frequency)

    if row[' Gene'] != ' intergenic' :
        if row[' Annotation'] == ' missense':
            aaBefore, aaAfter = utils_LTEE.geneModification(geneList, seq_ref, ref_codons, row[' Allele'][-1], row['Position'], strand)
            df_filtered.at[index, 'AA change'] = f'{aaBefore}->{aaAfter}'
            header = f'{row['Position']} {row[' Allele']} {row[' Annotation']} {aaBefore}->{aaAfter}'
        else:
            header = f'{row['Position']} {row[' Allele']} {row[' Annotation']}'
    else:
        header = f'{row['Position']} {row['Intergenic between']} {row[' Annotation']}'

    df_plot[header] = frequency_list


# Plotting with or without xkcd style and save it in output folder with a name containing population and gene of interest.
if xkcd == True:
    cm = plt.xkcd()
else:
    cm = nullcontext()

with cm :
    fig, axs = plt.subplots()
    
    df_plot.plot(ax = axs,
        marker = ".",
        linestyle = "-",
        title = f'{popName} - {gene}',
        figsize = (20, 5),
        fontsize = 16,
        xlim = (0, 60000),
        alpha = 0.5,
    )
    
    axs.legend(loc='best')
    
    axs.set_ylabel("Mutation frequency")
    axs.set_xlabel("Generation time")
    
    # Show or save the plot
    if plot_show == True :
        plt.show()
    
    else:
        fig.savefig(f'output/{pop}_{gene}.png')
        print(f'Plot saved in: {cwd}/output/{pop}_{gene}.png')

# Save datas in output folder as csv with a name containing population and gene of interest.
# The csv file contains all the information about mutation (position, allele, annotation, AA change if missense mutation) and the frequency for each time point.
if datas_save == True:
    if include_intergenic == True:
        new_col = df_columns[0:2] + ['Intergenic between'] + df_columns[2:4] + ['AA change'] + df_columns[4:]
        df_reordered = df_filtered.loc[:, new_col]
        df_reordered.to_csv(f'output/{pop}_{gene}_n_intergenic.csv', index=False, sep=',')
        print(f'Datas saved in: {cwd}/output/{pop}_{gene}_n_intergenic.csv')
    else:
        new_col = df_columns[0:4] + ['AA change'] + df_columns[4:]
        df_reordered = df_filtered.loc[:, ] 
        df_reordered.to_csv(f'output/{pop}_{gene}.csv', index=False, sep=',')
        print(f'Datas saved in: {cwd}/output/{pop}_{gene}.csv')