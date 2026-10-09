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

## 5. Ajustements sur SPARC (`examples/fit_sparc.py`)

163 galaxies de qualité 1-2 (sur 175), Upsilon_disque avec prior log-normale (0,5 ± 0,1 dex),
distances et inclinaisons non modélisées. Résultats du 9 octobre 2026.

### 5.1 Variantes de l'éq. (5)

La forme lue dans le manuscrit donne g < g_baryonique pour 23,7 % des points (r ≫ rc et g_bar ≫ a₀).
Deux variantes qui évitent ce défaut, avec w = x²/(1+x²) et x = r/rc :

| Variante | Rigidité | Couplage | g_eff | Propriétés |
|---|---|---|---|---|
| Manuscrit (lecture B) | R₀(1+x²)/(1+x² a₀/g) | racine carrée | g √((1+x² a₀/g)/(1+x²)) | g_eff < g pour r ≫ rc, g ≫ a₀ ; R_v peut dépasser R₀ |
| **V1** | R₀ / (1 + w a₀/g) | racine carrée | g √(1 + w a₀/g) | g_eff ≥ g ; R_v ≤ R₀ ; Newton si r ≪ rc ou g ≫ a₀ ; √(g a₀) si r ≫ rc, g ≪ a₀ |
| **V2** | R₀ / (1 + w √(a₀/g)) | G_eff = 1/R_v (cohérent RG) | g + w √(g a₀) | mêmes limites ; **le plateau vient de R_v sans racine carrée ajoutée à la main** |

### 5.2 Résultats

| Modèle | Paramètres libres par galaxie | χ²/dof médian | χ²/dof total | BIC | rms (dex) |
|---|---|---|---|---|---|
| Newton | Υ | 40,6 | 256,6 | 798 410 | 0,361 |
| MOND (a₀ fixé) | Υ | 3,8 | 14,2 | 45 502 | 0,282 |
| **MOND, a₀ libre (contrôle)** | Υ, a₀ | **2,0** | **5,3** | **18 297** | **0,114** |
| TCE manuscrit (rc libre) | Υ, rc | 3,5 | 20,4 | 62 514 | 0,126 |
| TCE V1 (rc libre) | Υ, rc | 2,8 | 12,4 | 39 045 | 0,117 |
| TCE V2 (rc libre) | Υ, rc | 2,3 | 7,8 | 25 556 | 0,117 |
| TCE manuscrit, rc = k Rdisk | Υ (+ k global = 1,26) | 4,9 | 30,5 | 96 072 | 0,165 |
| TCE V1, rc = k Rdisk | Υ (+ k global = 0,40) | 4,0 | 15,4 | 49 128 | 0,177 |
| TCE V2, rc = k Rdisk | Υ (+ k global = 0,40) | 4,5 | 15,4 | 49 135 | 0,188 |

### 5.3 Lecture prudente

- **Le défaut est corrigé.** V1 et V2 ne donnent jamais g < g_baryonique (contre 23,7 %).
- **V2 est la meilleure des formes TCE** : χ²/dof total 7,8 contre 14,2 pour MOND à a₀ fixé, et un BIC plus bas malgré un paramètre libre de plus par galaxie.
- **Mais le contrôle est décisif** : MOND avec a₀ libre par galaxie fait mieux que toutes les variantes TCE (χ²/dof total 5,3 ; BIC 18 297). L'avantage de TCE sur MOND vient donc de la liberté d'un paramètre par galaxie, pas de la forme de la loi. À flexibilité égale, ces données ne favorisent pas TCE.
- **rc = k·Rdisk avec un seul k** reste moins bon que MOND à a₀ fixé pour le manuscrit, et à peine meilleur pour V1 et V2 (χ²/dof total 15,4 contre 14,2) : rc n'est pas fixé par une échelle simple.
- Les variantes V1 et V2 sont des formes écrites après avoir vu le défaut : elles n'ont pas de dérivation, seulement des limites correctes.
- **Le problème des horloges (section 1) persiste.** Dans V1 et V2, R₀/R_v = 1 + w·a₀/g (ou 1 + w√(a₀/g)) atteint ~100 à faible accélération : si dt gouverne les horloges atomiques, l'écart reste d'ordre 1 là où la RG donne 10⁻⁶ à 10⁻⁸.
- Les χ² absolus ne sont pas interprétables sans les incertitudes de distance et d'inclinaison ; seules les comparaisons relatives le sont.
