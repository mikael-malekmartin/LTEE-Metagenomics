import pandas as pd
import numpy as np
import os

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import argparse
import textwrap

#=================
# ARGS
#=================
# Create an argument parser to get the population name, gene name and options for the plot and datas saving
parser = argparse.ArgumentParser(description=textwrap.dedent('''\
                                                Generate a circlar plot of mutations for a specific population along the LTEE.
                                                Default output folder for datas and plots : [APP DIR]/output/
                                                '''),
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

requiredNamed = parser.add_argument_group('required named arguments')
requiredNamed.add_argument("-p", "--population", help="Population name. Format : m/p1-6 (i.e. m1 or p6)", required=True)
requiredNamed.add_argument("-g", "--genes", help="A list of specific genes. Carefull with the case (i.e. topA spoT nadR)", nargs='+', default=[])
parser.add_argument("-c", "--clones", help="[Optionnal] Specifiy A or B if you want plot only mutations for a specific clone", nargs='+', default=['A', 'B'])
parser.add_argument("-a", "--annotations", action="store_true", help="[Optionnal] Add annotations for fixed mutated genes (carefull with the readability of the plot if too many annotations)")
parser.add_argument("-s", "--show_plot", action="store_true", help="[Optionnal] Show the plot, to edit before save it (specified destination)")

args = parser.parse_args()

config = vars(args)

pop = config.get('population').casefold()
gene = config.get('genes')
clone = config.get('clones')
clone = [c.upper() for c in clone]
annontations = config.get('annotations')
plot_show = config.get('show_plot')

os.makedirs('output', exist_ok=True) # Create folders if not exist for output (plot and datas)

goldenNumber = 1.61803398875

# 606 Génome size
genomeSize = 4629812

colours = ['#003f5c', '#2f4b7c', '#665191', '#a05195', '#d45087', '#f95d6a', '#ff7c43', '#ffa600'] # https://www.learnui.design/tools/data-color-picker.html

#=============
# PARAMETERS
#=============

pop_list = ['m1', 'm2', 'm3', 'm4', 'm5', 'm6', 
            'p1', 'p2', 'p3', 'p4', 'p5', 'p6']

if pop not in pop_list:
    raise ValueError(f"'{pop}' is not a correct population name. Try with the following format m/p1-6 (i.e. m1 or p6)")

# Convert the given argument 'pop' to Population name
minusOrPlus = pop[0]

if minusOrPlus == 'm' :
    minusOrPlus = '-'
elif minusOrPlus == 'p' :
    minusOrPlus = '+'

population = f'Ara{minusOrPlus}{pop[-1]}'


time_list = [606, 500, 1000, 2000, 5000, 10000, 20000, 50000]
mutation_type = ['SNP', 'MOB']   
#clone = ['A', 'B']
#annontations = True

oricStart = 3886105
tercStart = 1575774

specificGenes_bool = False

if len(gene) > 0 :
    specificGenes_bool = True
    specificGenes = gene

# Load data
try :
    df_raw = pd.read_csv('input/LTEE-Ecoli-data.csv', delimiter=",", index_col=0)
except FileNotFoundError as exc:
    raise FileNotFoundError(
        "Error: No file 'input/LTEE-Ecoli-data.csv' found. "
        "Check if the file is in the input folder."
    ) from exc

# Plotting

fig = plt.figure(figsize=(15,10),
                linewidth=10,
                )

rect = [0.1, 0.1, 0.8, 0.8]
offsetPlot = oricStart * 360 / genomeSize - 360

ax_polar = fig.add_axes(rect, polar=True, frameon=False)
ax_polar.set_theta_zero_location('N', offset=offsetPlot) # Offset en deg
ax_polar.set_theta_direction(-1)
ax_polar.set_title(f'Mutations overview for {population} [Clone {"&".join(clone)}] ({", ".join(mutation_type)})', x=1, y=1)
if specificGenes_bool != False:
    ax_polar.set_title(f'Mutations overview for {population} [Clone {"&".join(clone)}] ({", ".join(mutation_type)}) - {", ".join(specificGenes)}')

ax_polar.set_rgrids(range(len(time_list)), 
                    labels=time_list, 
                    angle=offsetPlot,
                    fontsize=14, fontweight='bold',
                    color='grey', verticalalignment='bottom', horizontalalignment='center')

circle_size = 2
initial_circle_height = 1
color_ring = 0

# Add sub divisions
for posT in range(0, genomeSize, 500000):
    x = posT
    ax_polar.annotate(x, xy = (x * circle_size * np.pi / genomeSize, len(time_list) - 0.7),
                textcoords='data',
                color='lightgrey',
                arrowprops=dict(arrowstyle="-", color="lightgrey" ),
                xytext=(x * circle_size * np.pi / genomeSize, len(time_list) + 0.75),
                horizontalalignment='center',
                verticalalignment='center')
    
for posTs in range(0, genomeSize, 100000):
    x2 = posTs
    ax_polar.annotate(text="", xy = (x2 * circle_size * np.pi / genomeSize, len(time_list) - 0.7),
                textcoords='data',
                color='lightgrey',
                arrowprops=dict(arrowstyle='-', color="lightgrey" ),
                xytext=(x2 * circle_size * np.pi / genomeSize, len(time_list) - 0.25),
                horizontalalignment='center',
                verticalalignment='center')

# Add OriC position
ax_polar.annotate(text="OriC", xy = (oricStart * circle_size * np.pi / genomeSize, 0),
    textcoords='data',
    color='green',
    arrowprops=dict(arrowstyle="-", ls="dashed", color="lightgrey" ),
    xytext=(oricStart * circle_size * np.pi / genomeSize, len(time_list) + 1),
    horizontalalignment='center',
    verticalalignment='center')

# Add TerC position
ax_polar.annotate(text="TerC", xy = (tercStart * circle_size * np.pi / genomeSize, 0),
    textcoords='data',
    color='green',
    arrowprops=dict(arrowstyle="-", ls="dashed", color="lightgrey" ),
    xytext=(tercStart * circle_size * np.pi / genomeSize, len(time_list) + 1),
    horizontalalignment='center',
    verticalalignment='center')

for time in time_list:
    ax_polar.barh(str(time), 
                genomeSize * circle_size * np.pi / genomeSize, 
                color=colours[color_ring], 
                height = goldenNumber * 0.4,
                alpha = goldenNumber * 0.15)
    color_ring += 1

j = 1
geneAnnotated = []
gene_list = {}
timeLen = 1
for time in time_list[1::] :
    df_time = df_raw[(df_raw['population'] == population) & (df_raw['time'] == time) & (df_raw['clone'].isin(clone)) & (df_raw['type'].isin(mutation_type))]
    if specificGenes_bool != False:
        df_time = df_time[df_time['gene_name'].isin(specificGenes)]

    mutationCount = 1

    for index, row in df_time.iterrows():
        if len(clone) > 1 and mutationCount > 1:
            mutationCount = 1
            continue
    
        gene = row['gene_name']
        cloneAB =  row['clone']
        generation = row['time']
        mutationPosition = row['start_position']

        if len(clone) > 1 and cloneAB == 'A':
            colorMut = 'red'
        elif len(clone) > 1 and cloneAB =='B':
            colorMut = 'green'
        else:
            colorMut = 'black'

        if len(clone) > 1 and len(df_time[df_time['start_position'] == mutationPosition]) > 1 :
            mutationCount = 2
            colorMut = 'black'

        mutation = f'{gene}-{mutationPosition}'
        
        if mutation in gene_list and generation != gene_list[mutation]['last'] :
            temp = gene_list[mutation]['value']
            temp += 2**timeLen
            gene_list.update({mutation : {'value' : temp,
                                        'gene': gene,
                                        'position': mutationPosition,
                                        'last' : generation}})
        elif mutation not in gene_list :
            gene_list.update({mutation : {'value' : 2**timeLen,
                                        'gene': gene,
                                        'position': mutationPosition,
                                        'last' : generation}})

        position = int(row['start_position'])
        ax_polar.barh(str(time),
                left = position * circle_size * np.pi / genomeSize,
                width = 1 * 0.005,
                color=colorMut,
                height = 1)
        
    timeLen += 1    
    j += 1

if (annontations == True) or (specificGenes_bool == True):
    df_annotation = pd.DataFrame.from_dict(gene_list, orient="index").sort_values(by='value', ascending=False)
    if (specificGenes_bool != True):
        df_annotation = df_annotation[(df_annotation['last'] == time_list[-1]) & (df_annotation['value'] > 2**(len(time_list)-1))]
    for index, row in df_annotation.iterrows():
        annotPositionMutation = row['position']
        annotMutatedGene = row['gene']

        ax_polar.annotate(annotMutatedGene, xy = (annotPositionMutation * circle_size * np.pi / genomeSize, 0),
            textcoords='data',
            color="black",
            #arrowprops=dict(arrowstyle= "-", ls= "dashed", color="grey"),
            xytext=(annotPositionMutation * circle_size * np.pi / genomeSize, len(time_list)),
            horizontalalignment='center',
            verticalalignment='center')
    
                

# Hide all grid elements
ax_polar.grid(False)
ax_polar.tick_params(axis='both', left=False, bottom=False, labelbottom=False, labelleft=True)

#Custom legend
if len(clone) > 1:
    legend_element = [Line2D([0], [0], color='black', lw=1),
                    Line2D([0], [0], color='red', lw=1),
                    Line2D([0], [0], color='green', lw=1)]
    
    plt.legend(legend_element, 
                ['Mutations in Clone A & B',
                'Mutations in Clone A',
                'Mutations in Clone B'],
                loc='upper right',
                bbox_to_anchor = (1.2, 1))


if plot_show == True:
    plt.show()
else :
    if specificGenes_bool != False:    
        fig.savefig(f'output/circular_plot_{population}-{"_".join(clone)}-{"_".join(specificGenes)}.png', bbox_inches='tight')
        print(f'Plot saved to output/circular_plot_{population}-{"_".join(clone)}-{"_".join(specificGenes)}.png')
    else:
        fig.savefig(f'output/circular_plot_{population}-{"_".join(clone)}.png', bbox_inches='tight')
        print(f'Plot saved to output/circular_plot_{population}-{"_".join(clone)}.png')
