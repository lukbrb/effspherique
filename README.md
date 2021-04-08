# Projet_Physique
Regroupe différents projets numériques de physique 

## Projet 1 ##

#### Le projet 1 est un travail sur l'effondrement à symétrie sphérique en cosmologie. ####

Y sont traités :

- Le temps d'effondrement, calculé à l'aide de Runge-Kutta 2 et Runge-Kutta 4
  - Fichier : Sol_nonlin.py
- Une vérification sur leur efficacité dans un régime linéaire (surdensité très inférieure à 1), avec une comparaison de la soltion calculée par scipy.integrate.odeint
  - Fichier : Sol_lin.py
- L'influence des différents paramètres initiaux sur le temps d'effondrement 
  - Fichiers : Sol_nonlin.py, Itération.py et les fichiers de résultats dans le dossier Résultats
- Le calcul du rayon de virialisation et de la surdensité viriel 
  - Fichier : Viriel.py
- Une animation permettant de mieux visualiser le projet
  - Fichiers : animation.mp4 et Graph_3D.py


#### A noter : Ceci était un premier projet et une découverte de Python, il est possible voire qu'il contienne certaines imperfections ####

### Mise à jour : un [notebook Jupyter](https://github.com/lukbrb/Projet_Physique/blob/master/Projet1-bis/ProjetNumerique.ipynb) a été crée, améliorant et éclaircissant les scripts déjà présents, et ajoutant des explications sur la démarche physique du projet
