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

Échantillon : 138 galaxies de qualité 1-2 ayant au moins 8 points valides (une galaxie qui a moins de
points que de paramètres libres n'a pas de degrés de liberté). Les chiffres de la version précédente
(163 galaxies) sont remplacés. Colonnes : χ²/dof médian, χ²/dof total, BIC, rms (dex).

**Sans nuisance (D et i fixés aux valeurs publiées)** — `python examples/fit_sparc.py`

| Modèle | Param. libres / galaxie | χ²/dof méd. | χ²/dof tot. | BIC | rms |
|---|---|---|---|---|---|
| Newton | Υ | 53,1 | 262,6 | 783 516 | 0,361 |
| MOND (a₀ fixé) | Υ | 4,1 | 14,6 | 44 574 | 0,286 |
| MOND, a₀ libre (contrôle) | Υ, a₀ | 2,4 | 5,4 | 17 678 | 0,114 |
| TCE manuscrit (rc libre) | Υ, rc | 3,7 | 20,9 | 61 689 | 0,126 |
| TCE V1 (rc libre) | Υ, rc | 3,1 | 12,7 | 38 272 | 0,117 |
| TCE V2 (rc libre) | Υ, rc | 2,5 | 8,0 | 24 867 | 0,118 |
| TCE V2, rc = k·Rdisk (k = 0,40) | Υ | 4,7 | 15,8 | 48 076 | 0,188 |

**Avec distance et inclinaison marginalisées** — `python examples/fit_sparc.py --marginalize`.
Pour chaque galaxie, la distance D' = D(1 + z·e_D/D) et l'inclinaison i' = i + z·e_i ont une prior
gaussienne (grille 7 × 7 sur ±2σ) ; les rayons deviennent r·D'/D, les vitesses baryoniques
V_bar·√(D'/D), les vitesses observées V_obs·sin i / sin i'. Le χ² contient les priors (maximum a posteriori).

| Modèle | Param. libres / galaxie | χ²/dof méd. | χ²/dof tot. | BIC | rms |
|---|---|---|---|---|---|
| Newton | Υ, D, i | 29,6 | 199,6 | 543 122 | 0,282 |
| MOND (a₀ fixé) | Υ, D, i | 3,0 | 7,0 | 22 189 | 0,278 |
| MOND, a₀ libre (contrôle) | Υ, D, i, a₀ | 2,6 | 5,6 | 18 760 | 0,111 |
| TCE manuscrit (rc libre) | Υ, D, i, rc | 3,0 | 11,6 | 34 162 | 0,107 |
| TCE V1 (rc libre) | Υ, D, i, rc | 2,4 | 6,9 | 22 152 | 0,100 |
| **TCE V2 (rc libre)** | Υ, D, i, rc | **2,0** | **4,5** | **16 037** | 0,100 |
| TCE V2, rc = k·Rdisk (k = 0,32) | Υ, D, i | 2,8 | 6,9 | 21 904 | 0,168 |

### 5.3 Lecture prudente

- **Le défaut est corrigé.** V1 et V2 ne donnent jamais g < g_baryonique (contre 23,5 % pour le manuscrit).
- **Modéliser D et i change beaucoup.** Le χ²/dof total de MOND passe de 14,6 à 7,0 : une grande partie de l'écart apparent venait des incertitudes de distance et d'inclinaison. Les χ² absolus restent loin de 1 (autres systématiques non modélisées : asymétries, supports de pression, etc.).
- **Le classement dépend du traitement.** Sans nuisance, MOND à a₀ libre bat TCE V2 (5,4 contre 8,0). Avec nuisance, c'est TCE V2 qui bat MOND à a₀ libre (4,5 contre 5,6 ; BIC 16 037 contre 18 760). La conclusion « TCE n'apporte rien à flexibilité égale » du premier ajustement n'est donc pas robuste.
- **Mais le contrôle est imparfait.** Les deux familles à un paramètre ne sont pas équivalentes : a₀ libre ne fait que décaler une transition de forme fixe, rc libre déplace la transition de V2 en rayon. Un contrôle équitable utiliserait une famille MOND à deux échelles (a₀ et un rayon de transition). Je ne l'ai pas fait.
- **rc = k·Rdisk avec un seul k** reste équivalent à MOND à a₀ fixé (χ²/dof total 6,9 contre 7,0) : l'avantage de TCE vient entièrement du rc propre à chaque galaxie.
- **Le problème des horloges (section 1) persiste.** Dans V1 et V2, R₀/R_v atteint ~100 à faible accélération : si dt gouverne les horloges atomiques, l'écart reste d'ordre 1 là où la RG donne 10⁻⁶ à 10⁻⁸.

## 6. Une loi pour rc à partir de la densité de surface baryonique ? (`examples/fit_sparc_laws.py`)

Lois testées sur TCE V2, distances et inclinaisons marginalisées, 138 galaxies, 3 118 points. Les paramètres
globaux sont ajustés sur une moitié des galaxies et évalués sur l'autre (validation croisée à 2 plis,
moyenne de 20 tirages). Σ_b = Υ*·SB_disque, Σ† = a₀/G ≈ 861 M☉/pc², M_bar = Υ*·L[3.6] + 1,33·M_HI.

| Modèle | Paramètres globaux | χ² total | χ² hors échantillon |
|---|---|---|---|
| MOND (a₀ fixé) | 0 | 18 858 | 18 858 |
| MOND, a₀ libre par galaxie | 0 (+1/galaxie) | 14 319 | — |
| TCE V2, rc libre par galaxie | 0 (+1/galaxie) | 11 596 | — |
| L1 : rc = k·Rdisk | 1 (k = 0,32) | 18 565 | 19 239 |
| L2 : rc = k·√(G M_bar / a₀) | 1 (k = 0,32) | 18 900 | 19 575 |
| L3 : rc = k·Rdisk·(Σ_b/Σ†)^α | 2 (k = 0,32 ; α = 0,00) | 18 565 | 19 801 |

**Résultat : aucune loi ne bat MOND hors échantillon.** L3 retombe sur L1 (α = 0, Δχ² = 0) : la densité
de surface n'apporte aucune information au-delà de Rdisk. Les paramètres optimaux (k ≈ 0,3) rendent
rc petit devant les rayons mesurés, et V2 tend alors vers g + √(g·a₀) sur toute la courbe, c'est-à-dire vers
une interpolation de type MOND. Tout le gain de rc libre (18 858 → 11 596) vient de la dispersion de rc d'une galaxie à l'autre,
que ni Rdisk, ni M_bar, ni Σ_b ne prédisent.

Corrélations de log rc libre (95 galaxies, rc hors bords de grille) : M_bar +0,63, Vflat +0,63 (circulaire : Vflat
vient de la courbe elle-même), Rdisk +0,52, Reff +0,50, Σ_b +0,42, fraction de gaz −0,45. Régression
log rc ~ Rdisk + Σ_b + M_bar : R² = 0,40, dispersion résiduelle 0,51 dex pour un intervalle de rc de 2,8 dex.
rc est donc en partie corrélé à la masse et à la taille, mais pas au point de définir une loi prédictive.

Conséquence pour la théorie : en l'état, rc est un paramètre libre par galaxie, pas une prédiction. Pour que
TCE soit prédictive, il faudrait une dérivation qui lie rc à des quantités mesurables indépendamment de la
courbe de rotation, ou un test d'une famille MOND à deux échelles pour savoir si rc fait vraiment mieux.

## 7. Contrôle équitable : TCE contre des familles MOND à un paramètre libre (`examples/fit_sparc_controls.py`)

Question laissée ouverte en 5.3 : l'avantage de rc libre vient-il de TCE ou de la présence d'un second paramètre ?
Chaque modèle a les mêmes nuisances et **un seul paramètre libre de plus par galaxie** (138 galaxies) :
`mond_a0` (a₀ libre), `mond_rs` (a₀ fixé, rayon de transition libre : g = g_bar + w·(g_MOND − g_bar), même structure que
TCE V2), `mond_n` (a₀ fixé, indice d'interpolation n libre), et TCE avec rc libre.

| Modèle (1 paramètre libre / galaxie) | χ² total, sans nuisance | BIC | χ² total, D et i marginalisés | BIC |
|---|---|---|---|---|
| MOND (a₀ fixé, 0 paramètre) | 43 464 | 44 574 | 18 858 | 22 189 |
| MOND, a₀ libre | **15 457** | **17 678** | 14 319 | 18 760 |
| MOND, rayon de transition libre (`mond_rs`) | 26 372 | 28 592 | 13 644 | 18 085 |
| MOND, indice d'interpolation libre (`mond_n`) | 29 461 | 31 682 | 14 937 | 19 378 |
| TCE du manuscrit, rc libre | 59 469 | 61 689 | 29 721 | 34 162 |
| TCE V1, rc libre | 36 051 | 38 272 | 17 711 | 22 152 |
| TCE V2, rc libre | 22 646 | 24 867 | **11 596** | **16 037** |

Test de signe par galaxie (même nombre de paramètres), D et i marginalisés :

| A contre B | A meilleur sur | p | médiane Δχ² |
|---|---|---|---|
| TCE V2 contre `mond_rs` | 85 / 138 | 0,008 | −0,8 |
| TCE V2 contre `mond_a0` | 95 / 138 | 1×10⁻⁵ | −3,9 |
| TCE V2 contre `mond_n` | 100 / 138 | 1×10⁻⁷ | −3,2 |
| TCE V1 contre `mond_rs` | 39 / 138 | 3×10⁻⁷ | +1,5 |

Sans nuisance, les mêmes tests donnent TCE V2 à égalité avec `mond_rs` (73/138, p = 0,55) et en retrait de `mond_a0`.

### Lecture prudente

- **Une partie de l'avantage de rc libre est générique** : MOND avec n'importe quel second paramètre passe de χ²/dof 7,0 à 5,3-5,8 une fois D et i marginalisés. Seul un contrôle à nombre de paramètres égal permet de le voir.
- **TCE V2 reste devant, de peu.** Δχ² = 2 048 sur `mond_rs`, mais seulement 0,8 de médiane par galaxie ; l'effet est statistiquement significatif (p = 0,008) quand D et i sont marginalisés, mais il disparaît (p = 0,55) quand ils ne le sont pas. Il est donc fragile devant les systématiques de modélisation.
- **D'où vient la différence de forme ?** V2 écrit g = g_bar + w·√(g_bar·a₀), qui revient vers Newton comme √(a₀/g_bar), plus lentement que l'interpolation « simple » de MOND (a₀/g_bar). Les données semblent légèrement préférer ce retour lent, mais `mond_n` (indice libre) ne le reproduit pas, ce qui mérite une étude séparée.
- **V1 est moins bon que `mond_rs`** (39/138 galaxies) : l'avantage de V2 ne vient pas de la structure en « portail », mais de la forme précise de la transition.
- **V2 a été écrit après avoir vu le défaut de la forme (5)** ; choisir la forme parmi plusieurs candidats puis la comparer sur les mêmes données gonfle l'avantage apparent. Il faudrait la tester sur des données indépendantes (autre échantillon que SPARC).
- Les remarques précédentes tiennent toujours : rc reste un paramètre libre par galaxie (section 6) et le problème des horloges (section 1) n'est pas résolu.

## 8. Test hors échantillon sur THINGS / LITTLE THINGS

### 8.1 Ce qui a pu être fait, et ce qui ne l'a pas pu

| Élément | Statut |
|---|---|
| Courbes de rotation **baryoniques radiales** de THINGS (de Blok+2008 : Vgas, Vstar) | **Non récupérées** : la page de données de l'enquête THINGS ne contient que des cubes et cartes FITS bruts, et la table de de Blok+2008 n'est pas sur le CDS. Un ajustement radial de TCE V2 sur THINGS exige ces courbes |
| LITTLE THINGS (Oh+2015, CDS J/AJ/149/180) | Récupéré : 26 naines HI, mais le CDS ne fournit que des courbes de rotation totales normalisées et des masses intégrées (gaz, étoiles), pas de profils baryoniques radiaux |
| Test de la relation de Tully–Fisher baryonique (éq. 7) sur LITTLE THINGS | **Fait** (`examples/test_littlethings_btfr.py`) |
| Chaîne d'ajustement radial prête pour un échantillon externe | **Prête** (`sparc.load_rotmod_directory` et `examples/fit_external.py`), testée sur des données synthétiques seulement |

Trois galaxies de LITTLE THINGS sont déjà dans SPARC (DDO 154, DDO 168, NGC 2366) et sont exclues ; 4 autres n'ont pas de masse stellaire. Reste 19 galaxies indépendantes.

### 8.2 Résultats de la relation de Tully–Fisher baryonique, v⁴ = G·M_bar·a₀ (a₀ = 1,2×10⁻¹⁰ m/s²)

| Échantillon | N | ⟨log(v_obs/v_pred)⟩ | a₀ ajusté (m/s²) | Pente de la BTFR [16-84 %] |
|---|---|---|---|---|
| SPARC (Vflat, qualité 1) | 87 | +0,037 ± 0,007 | 1,7×10⁻¹⁰ | 3,56 [3,45–3,68] |
| LITTLE THINGS, V(Rmax) | 19 | +0,009 ± 0,043 | 1,3×10⁻¹⁰ | 1,68 [0,98–2,18] |
| LITTLE THINGS, Viso(Rmax) | 19 | +0,021 ± 0,037 | 1,5×10⁻¹⁰ | 1,98 [1,42–2,43] |
| LITTLE THINGS, gaz dominant | 11 | −0,059 ± 0,040 | 7×10⁻¹¹ | 2,65 [1,50–3,30] |

Lecture prudente :
- **La normalisation est compatible** : sur un échantillon indépendant de SPARC, la vitesse prédite par l'éq. (7) avec a₀ = 1,2×10⁻¹⁰ m/s² est correcte à 0,01 ± 0,04 dex près (a₀ ajusté 1,3×10⁻¹⁰), donc le même a₀ vaut pour les naines. C'est le résultat que MOND prédit aussi : ce test ne distingue pas TCE de MOND.
- **La pente n'est pas reproduite** : 1,7 [1,0–2,2] au lieu de 4, et même la pente de SPARC (3,56) est inférieure à 4 avec cette recette. L'échantillon LITTLE THINGS couvre seulement 1 dex en vitesse (12–126 km/s), la masse varie de 0,5 dex à 1 dex pour une même vitesse (dispersion 0,5 dex), et la pente dépend fortement d'une galaxie : DDO 210 (V = 12 km/s) la fait varier de 1,05 à 1,88 quand on la retire.
- **V(Rmax) n'est pas une vitesse plate** pour beaucoup de naines (la courbe monte encore), ce qui biaise la pente vers le bas. Un test propre de la pente demande une vitesse plate bien définie.
- **Ce test ne compare pas les formes de rotation** (V1, V2, MOND) : il ne valide donc pas la forme V2, écrite après avoir vu SPARC. C'est le test radial sur des courbes baryoniques indépendantes qui le ferait.

### 8.3 Pour faire le test radial sur THINGS

Il faut fournir les courbes baryoniques radiales (gaz et étoiles, rapport masse/luminosité à 1) de THINGS au format SPARC `*_rotmod.dat`
(colonnes `Rad Vobs errV Vgas Vdisk Vbul SBdisk SBbul`, en-tête `# Distance = X Mpc`), puis lancer :

```
python examples/fit_external.py --dir data/things_rotmod
```

Le script exclut les galaxies déjà dans SPARC et compare Newton, MOND (a₀ fixé, a₀ libre, rayon de transition libre, indice libre) et TCE V1/V2 à un paramètre libre par galaxie. Seul un résultat où V2 bat `mond_rs` sur ces données indépendantes, avec le même nombre de paramètres, validerait l'avantage mesuré sur SPARC (section 7).
