# Projets de Physique

## 1 - Effondrement à symétrie sphérique en cosmologie

Dans ce projet nous allons déterminer le temps d’effondrement d’une surdensité sphérique
dans un univers en expansion. Dans le cadre du modèle standard de la cosmologie, ce sont les
surdensités primordiales qui, en s’effondrant, conduisent à la formation des halos de matière noire
au sein desquelles le gaz tombe pour former des galaxies. En connaissant ce temps d’effondrement
on peut donc en principe déterminer le rythme de formation des galaxies et amas de galaxies
ainsi que leur abondance dans notre univers.
Nous nous intéressons à la résolution numérique du problème, à savoir les multiples types de méthodes d’intégration et leur pertinence. Nous étudions également l’influence des différents paramètres physique ou numérique sur la solution finale. Finalement, nous faisons un retour sur certaines notions de
physique (la surdensité viriel) qui permettent de rendre notre étude plus réaliste.

Les solutions ont initiallement été implémentées en Python, puis adaptées en Fortran pour gagner en efficacité, ainsi que pour comparer les temps d'éxécution. Le code en Python se trouve dans le dossier `effspher-py` tandis que le code en Fortran se trouve dans le dossier `effspher-f`.

**Prérequis et installation :** Pour le code en Python, seules les bibliothèques numpy et matplotlib sont utilsées, fournies directement avec une installation Anaconda/Miniconda. Pour le Fortan, l'utilisation du gestionnaire de paquet et de compilation [fpm](https://fpm.fortran-lang.org/fr/index.html) est recommandée.
