---
title: "Point process discrimination according to repulsion"
author: "H. Adrat and L. Decreusefond"
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---
# Point process discrimination according to repulsion

## Abstract

In numerous applications, cloud of points do seem to exhibit {repulsion} in the intuitive sense that there is no local cluster as in a Poisson process. Motivated by data coming from cellular networks, we devise a classification algorithm based on the form of the Voronoi cells. We show that, for data  representing the locations of antennas in a cellular networks, we can retrieve some repulsion between antennas, a feature that  was expected for engineering reasons.
 
# Introduction

In the performance analysis of cellular systems, the locations of antennas (or base stations) play a major role (see {cite}`BaccelliStochasticGeometryWireless2008`). It is usually admitted that they can be modeled by a Poisson process. But the data which can be gathered from the Web site of the French National Agency of Radio Frequencies, Cartoradio, see {cite}`ANFR`, tend to prove that this may not be the case. More precisely, if we look at the global picture of all antennas in Paris, we see features reminiscent of a Poisson process (local clusters for instance), see Figure~\ref{fig:antennas}. However, if we look closer and finer, by specifying a region and a frequency band, we see that the antennas locations do seem to exhibit some repulsion (see Figure~\ref{fig:antennas}, right picture). 

````{prf:definition}
Es sei $U\subset \C$ offen und seien $\gamma_0,\gamma_1:[a,b]\to U$ zwei Wege. Falls eine Zerlegung $(a=t_0,\ldots,b=t_N)$ existiert mit Kreisscheiben $B^j$, welche zulässig für beide Wege ist, d.h.,

```{math}
\gamma_0([t_j, t_{j+1}])\cup \gamma_1([t_j,t_{j+1}]) \subset B^j
```

dann nennen wir $\gamma_0$ und $\gamma_1$ **benachbart**.
````



```{bibliography}
:style: unsrt
```
