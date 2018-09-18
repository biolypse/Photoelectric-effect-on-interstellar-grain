Pour générer un grain:
- ouvrir un terminal
- se placer dans le répertoire contenant GenGrain.py
- taper la commande : 
   > time python GenGrain.py

Rq: 
- La commande "time" permettra de mesurer le temps d'exécution du script
- Il est possible d'ouvrir GenGrain.py avec un éditeur de texte et de changer les paramètres. Pensez à conserver le script d'origine !
- Le script GRF_routines.py doit être présent dans le même répertoire que GenGrain.py car il est appelé par GenGrain.py. Il est inutile (déconseillé) d'ouvrir GRF_routines.py.
- Il est possible de choisir de visualiser ou non les grains au moment où ils sont produits avec la variable "doplot" dans GenGrain.py.
- Le grain est écrit dans un fichier ASCII nommé "Grain_***.txt", où les paramètres principaux sont notés.
