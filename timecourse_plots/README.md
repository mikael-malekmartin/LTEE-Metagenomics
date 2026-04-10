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
	python3 timecourse_plots.py -p POPULATION -g GENE [-s] [-i] [-d]

	options:
  		-h, --help            show this help message and exit
  		-i, --include_intergenic	[Optionnal] Include intergenic mutation before or after gene of interest
  		-s, --show_plot       		[Optionnal] Show the plot, to edit before save it (specified destination)
  		-d, --datas_save      		[Optionnal] Save datas as csv coma separate in default output folder

	required named arguments:
  		-p, --population POPULATION     Population name. Format : m/p1-6 (i.e. m1 or p6)
  		-g, --gene GENE       		The name of a gene (i.e. topA)
