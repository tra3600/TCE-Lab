# Chrono-Energetics : confrontation avec la RG et parallèles avec les théories émergentes

Ce document accompagne `tce/confrontation.py` et `examples/chrono_energetics_demo.py`
(étapes 6 à 9). Les résultats numériques viennent du code ; les références
viennent d'une recherche web faite pendant la session (sources secondaires :
résumés et pages d'accueil d'articles, pas les articles complets sauf mention).
Les mentions « connaissance générale » ne sont pas vérifiées dans cette session.

## 1. Confrontation avec le modèle relativiste

| Test | Résultat dans le code | Lecture |
|---|---|---|
| Trou noir de Schwarzschild | `R_v = R_0/√(1−r_s/r)` redonne dτ = √(1−r_s/r)·dt | Reformulation de la RG par construction, pas une prédiction |
| Horloges en galaxie (éq. 1 + éq. 5) | Écart au taux 1 : −10 % à 1 kpc, +670 % à 30 kpc, ×86 à 100 kpc ; RG : ~10⁻⁶ à 10⁻⁸ | Si dt gouverne les horloges atomiques, l'écart est de 5 à 8 ordres de grandeur trop grand. Il faut préciser ce que dt gouverne |
| Couplage g_eff | `G_eff = 1/R_v` (lecture RG) donne une courbe qui **monte** (608 km/s à 100 kpc) ; la racine carrée du manuscrit donne le plateau | La racine carrée de l'éq. (6) n'est pas dérivée de R_0 = 1/κ ; c'est un ajout à justifier |
| Courbes de rotation | TCE (racine), MOND et RG+NFW sont proches à 10–100 kpc ; seul TCE (et MOND) donne un plateau exact et v⁴ = G M a₀ | Le modèle ne se distingue de MOND que par la dépendance en `rc` |
| Lentilles gravitationnelles | Non calculable : le manuscrit n'a pas de métrique, donc pas de déflexion de la lumière | Condition nécessaire pour tout concurrent sérieux de la matière noire |
| Échelle a₀ | cH₀/2π = 1,08×10⁻¹⁰, cH₀/6 = 1,13×10⁻¹⁰, mesuré ~1,2×10⁻¹⁰ m/s² (H₀ = 70 km/s/Mpc) | Coïncidence déjà notée en MOND ; TCE la traite comme paramètre libre |

## 2. Parallèles avec les théories émergentes et les travaux récents

| Cadre | Parallèle avec TCE | Différence | Statut (d'après les sources) |
|---|---|---|---|
| **Jacobson (1995)**, [entanglement equilibrium](https://ar5iv.arxiv.org/html/1810.12236) | R₀ = 1/κ = c⁴/(8πG) et R₀·c = ħc²η/2π avec η = 1/(4 l_P²) (étape 9 du code, test passant) | Jacobson dérive la RG et ne donne pas de loi dt = E/R | Cadre établi ; η est une entrée |
| **Sakharov (1967)**, [gravité induite](https://en.wikipedia.org/wiki/Induced_gravity) | « Rigidité du vide » = coefficient d'Einstein–Hilbert induit par les fluctuations quantiques | Sakharov part de champs quantiques, TCE postule R_v | Idée ancienne ; cadrage récent : Sakharov explique l'origine du coefficient, Jacobson pourquoi c'est une équation d'état |
| **Verlinde**, gravité émergente, [Brouwer et al. 2017](https://www.arxiv-vanity.com/papers/1612.03034/) | Même ambition : phénomène MOND comme effet entropique ; a₀ ~ cH₀ | Verlinde propose une dérivation (entropie en loi de volume) ; TCE postule le profil | Premier test par lentille faible : encourageant mais conclusion jugée prématurée par les auteurs |
| Critique de Verlinde : [Milgrom & Sanders](https://ar5iv.labs.arxiv.org/html/1612.09582) | Montre les pièges à éviter pour TCE | — | Formule jugée en conflit possible avec contraintes du système solaire et avec le cœur des amas |
| **Bianconi, « Gravity from entropy »**, PRD 111, 066001 (2025), [résumé QMUL](https://www.seresearch.qmul.ac.uk/news/4918/gravity-from-entropy-a-radical-new-approach-to-unifying-quantum-mechanics-and-general-relativity) | Gravité obtenue comme entropie relative quantique entre métriques ; loi d'aire de Schwarzschild | Pas de temps émergent ; TCE n'a pas d'action | Publié ; constante cosmologique émergente rapportée par presse (non vérifiée ici) |
| **Page–Wootters, Moreva et al.**, [expérience](https://arxiv.org/pdf/1310.4691) | Le temps comme corrélation (horloge-photon) vu par un observateur interne | TCE lie le temps à l'énergie dissipée, pas à l'intrication | Illustration expérimentale (2013-14), suivie d'un travail en 2026 sur l'accès partiel de l'observateur |
| **Temps thermique, Connes–Rovelli** (connaissance générale) | Temps et flot thermodynamique : proche de dt = dE/R si R joue le rôle d'une température | Pas d'état thermique explicite dans TCE | Hypothèse connue, non vérifiée ici |
| **MOND relativiste AeST**, [Skordis & Złośnik](https://ar5iv.arxiv.org/html/2007.00082) | Seul modèle de type MOND avec une action covariante, CMB ajusté, c_GW = c | AeST a une action et des champs ; TCE n'en a pas | [Panorama APS](https://physics.aps.org/articles/v14/143) ; la question de l'amas de la Balle n'est pas tranchée dans les sources trouvées ; critiques (« baroque ») |
| **Binaires larges**, [synthèse](https://www.universetoday.com/165136/the-debate-continues-do-wide-binaries-prove-or-disprove-mond) | Test direct du régime g ≪ a₀ à l'échelle stellaire : TCE (racine) et MOND y prédisent un excès de gravité | — | Débat ouvert : Chae et al. annoncent un facteur ~1,4-1,6 ; Banik et al. concluent à Newton ; la sensibilité à la méthode est soulignée (2026) |
| **Janus (Petit)**, voir PDF §7 | Deux secteurs, masse négative et temps inversé | Bimétrique ; critiqué | Minoritaire |

## 3. Ce qui manque à TCE pour entrer dans cette conversation

1. **Une action ou une dérivation** de R_v(r) à partir d'une loi microscopique (Jacobson + entropie en volume de Verlinde, ou Sakharov), au lieu du profil (5).
2. **Un couplage défini** entre R_v et la matière : qu'est-ce que dt gouverne (horloges atomiques, désintégrations, rayonnement) ? Le test d'horloges (étape 6) montre que la réponse contraint tout le reste.
3. **Une dérivation de la racine carrée** de l'éq. (6) (étape 7) : avec G_eff = 1/R_v on ne trouve pas de plateau.
4. **Un secteur relativiste** (métrique, lentilles, ondes gravitationnelles, CMB) comme dans AeST.
5. **Une confrontation aux données** : courbes de rotation SPARC (rc ajustable galaxie par galaxie), relation de Tully–Fisher baryonique, binaires larges, amas de la Balle avec un terme de transport (voir l'extension `advect_with_inertia`).
6. **L'amas de la Balle** : l'éq. (8) seule fait retarder la lentille sur le gaz (test `test_eq8_alone_lags_the_gas`).

## 4. Prochaines étapes proposées dans le dépôt

- Ajuster `rc` sur des courbes SPARC (téléchargement hors de cette session).
- Remplacer le profil (5) par la forme dérivée d'une entropie en loi de volume et comparer.
- Ajouter un test statistique sur les binaires larges (modèle de boost g_eff/g_N).

## 5. Premier ajustement sur SPARC (`examples/fit_sparc.py`)

163 galaxies de qualité 1-2 (sur 175), Upsilon_disque avec prior log-normale (0,5 ± 0,1 dex),
distances et inclinaisons non modélisées. Résultats du 9 octobre 2026 :

| Modèle | χ²/dof médian | χ²/dof total | BIC | rms (dex) |
|---|---|---|---|---|
| Newton (baryons) | 40,6 | 256,6 | 798 410 | 0,361 |
| MOND (a₀ fixé) | 3,8 | 14,2 | 45 502 | 0,282 |
| TCE (rc libre par galaxie) | 3,5 | 20,4 | 62 514 | 0,126 |
| TCE (rc = k·Rdisk, k = 1,26) | 4,9 | 30,5 | 96 072 | 0,165 |

Lecture prudente :
- TCE avec rc libre a le plus petit écart médian et le plus petit rms, mais il a un paramètre de plus par galaxie ; son χ² total et son BIC sont moins bons que ceux de MOND.
- Pour 23,7 % des points, TCE donne g < g_baryonique (jamais pour MOND) : c'est la limite connue de la forme (5) pour r ≫ rc et g_bar ≫ a₀ (gravité affaiblie). Elle rend le modèle non physique dans ces zones.
- rc libre est corrélé avec Rdisk (corrélation log-log 0,44), mais 48 galaxies sur 163 butent sur le bord de la grille de rc ; rc = k·Rdisk, avec un seul k, est nettement moins bon.
- Les χ² absolus ne sont pas interprétables sans les incertitudes de distance et d'inclinaison ; seules les comparaisons relatives le sont.
