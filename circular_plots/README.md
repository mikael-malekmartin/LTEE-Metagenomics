===============
Installation :
===============

1 - ouvrir un terminal dans le dossier du module

2 - exécuter la commande suivante :
	pip install -r requirements.txt

=========
Usage :
=========

1 - ouvrir un terminal dans le dossier du module

2 - exécuter la commande suivante :
	python3 -p POPULATION [-g GENES [GENES ...]] [-c CLONES [CLONES ...]] [-a] [-s]

	Generate a circlar plot of mutations for a specific population along the LTEE.
	Default output folder for datas and plots : [APP DIR]/output/

	options:
	-h, --help            	show this help message and exit
	-c, --clones CLONES [CLONES ...]
							[Optionnal] Specifiy A or B if you want plot only mutations for a specific clone
	-a, --annotations     	[Optionnal] Add annotations for fixed mutated genes (carefull with the readability of the plot if too many annotations)
	-s, --show_plot       	[Optionnal] Show the plot, to edit before save it (specified destination)

	required named arguments:
	-p, --population POPULATION
							Population name. Format : m/p1-6 (i.e. m1 or p6)
	-g, --genes GENES [GENES ...]
                     		A list of specific genes. Carefull with the case (i.e. topA spoT nadR)
