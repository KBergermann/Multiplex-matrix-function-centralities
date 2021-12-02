# Multiplex-matrix-function-centralities

This repository contains Matlab codes reproducing the results and grahpics from the

**Paper:** [1] Bergermann, K. & Stoll, M. (2021) Fast computation of matrix function-based centrality measures for layer-coupled multiplex networks. https://arxiv.org/abs/2104.14368v3

**Abstract:** Centrality measures identify and rank the most influential entities of complex networks. In this paper, we generalize matrix function-based centrality measures, which have been studied extensively for single-layer and temporal networks in recent years to layer-coupled multiplex networks. The layers of these networks can reflect different relationships and interactions between entities or changing interactions over time. We use the supra-adjacency matrix as network representation, which has already been used to generalize eigenvector centrality to temporal and multiplex networks. With a suitable choice of edge weights, the definition of single-layer matrix function-based centrality measures in terms of walks on networks carries over naturally to the multilayer case. In contrast to other walk-based centralities, matrix function-based centralities are parameterized measures, which have been shown to interpolate between (local) degree and (global) eigenvector centrality in the single-layer case. As the explicit evaluation of the involved matrix function expressions becomes infeasible for medium to large-scale networks, we present highly efficient approximation techniques from numerical linear algebra, which rely on Krylov subspace methods, Gauss quadrature, and stochastic trace estimation. We present extensive numerical studies on synthetic and real-world multiplex transportation, communication, and collaboration networks. The comparison with established multilayer centrality measures shows that our framework produces meaningful rankings of nodes, layers, and node-layer pairs. Furthermore, our experiments corroborate the linear computational complexity of the employed numerical methods in terms of the network size that is theoretically indicated under the assumption of sparsity in the supra-adjacency matrix. This excellent scalability allows the efficient treatment of large-scale networks with the number of node-layer pairs of order 10^7 or higher.


This repository contains:

**License:**
 - LICENSE: GNU General Public License v2.0

**Directories:**
 - Example_6_1_Synthetic_undirected_network: contains codes to reproduce [1, Fig. 5].
 - Example_6_2_Synthetic_temporal_network: contains codes to reproduce [1, Fig. 6].
 - Example_6_3_Numerical_approximation_error: contains codes to reproduce [1, Fig. 7].
 - Example_6_3_Numerical_approximation_error/data: contains Scotland Yard and the department 3 subset of the Email-EU [3] intra-layer adjacency matrices.
 - Example_6_4_European_airlines_network: contains codes to reproduce [1, Tab. 2-3] and [1, Fig. 8-9].
 - Example_6_4_European_airlines_network/data: contains intra-layer adjacency matrices as well as the node and layer ID lists of the largest connected cluster of the European airlines data set [2,3].
 - Example_6_5_Temporal_IMDb_Sci-Fi_network: contains codes to reproduce [1, Tab. 4] and [1, Fig. 10].
 - Example_6_5_Temporal_IMDb_Sci-Fi_network/data: contains intra-layer adjacency matrices as well as the node and layer ID lists of the temporal IMDb Sci-Fi network built by the authors from the full IMDb data set downloaded from https://datasets.imdbws.com.
 - Example_6_5_Temporal_IMDb_Sci-Fi_network/results: empty directory into which diagonal estimation results are saved for the computation of more accurate Gauss quadrature approximations for the top-ranked node-layer pairs.
 - subroutines: contains the subroutines required to reproduce the paper results (see 'Functions' below).
 - subroutines/funm_kryl: codes of funm_kryl Matlab toolbox for the approximation of f(A)b using restarted Krylov subspace methods (see [4] or the file subroutines/funm_kryl/index.htm).
 
 **Data:**
 - Example_6_3_Numerical_approximation_error/data/email_adjacencies_15d.mat: intra-layer adjacency matrices of a temporal network constructed from the department 3 subset of the Email-EU data set [5] in which each layer contains the number of exchanged e-mails between the employers of a European research institute as directed edges over a 15 day time interval.
 - Example_6_3_Numerical_approximation_error/data/scotlandyard_adjacency.mat: unweighted intra-layer adjacency matrices of a transportation network built from the Scotland Yard board game by the authors. Nodes represent stops in the city of London and layers represent different modes of transportation (boat, underground, bus and taxi).
 - Example_6_3_Numerical_approximation_error/data/scotlandyard_adjacency_weighted.mat: weighted intra-layer adjacency matrices of a transportation network built from the Scotland Yard board game by the authors. Nodes represent stops in the city of London and layers represent different modes of transportation (boat, underground, bus and taxi). Weights of the boat, underground, and bus layer are constructed from the minimal number of taxi rides required to travel the same distance.
 - Example_6_4_European_airlines_network/data/multiplex_airlines_airport_names_coords.mat: contains a list of the full airport names represented by the physical nodes in the European airlines network.
 - Example_6_4_European_airlines_network/data/multiplex_airlines_GC.mat: contains unweighted intra-layer adjacency marices of the largest connected cluster of the European airlines data set [2,3] as well as airport tokens as node IDs.
 - Example_6_4_European_airlines_network/data/multiplex_airlines_layer_ids.mat: contains a list of the airline names represented by the layers of the European airlines network.
 - Example_6_5_Temporal_IMDb_Sci-Fi_network/data/names_sci_fi.csv: contains a table of all principal cast/crew members of Sci-Fi productions included in the IMDb data set (including name IDs, real names, and additional information).
 - Example_6_5_Temporal_IMDb_Sci-Fi_network/data/sci_fi_seasons_intra_layer_adjacency.mat: intra-layer adjacency matrices of the temporal IMDb Sci-Fi collaboration network (files name.basics, title.basics, title.episode, and title.principals of the full IMDb data set downloaded from https://datasets.imdbws.com on 27.10.2021 14:00 CEST) as well as node and layer ID tables.
  

**Scripts:**
 - Example_6_1_Synthetic_undirected_network/small_synthetic_undiected_multiplex_example.m: creates a small synthetic multiplex network as well as its aggregated single-layer version and computes degree as well as matrix function-based centrality measures for both networks.
 - Example_6_2_Synthetic_temporal_network/synthetic_temporal_network_example.m: creates a small synthetic temporal network according to [6, Sec. 5.2] and compares dynamic communicabilities with marginal node Katz centralities for forward and backward evolving time.
 - Example_6_3_Numerical_approximation_error/email_KC.m: compares the approximation of Katz centrality obtained with Krylov subspace methods to the 'exact' quantity computed with the backslash for the temporal E-Mail EU network.
 - Example_6_3_Numerical_approximation_error/email_SC_bipartite.m: compares the upper bound on subgraph centrality computed by Gauss--Radau quadrature applied to the symmetric bipartite network representation to the 'exact' quantity computed with expm for the temporal E-Mail EU network.
 - Example_6_3_Numerical_approximation_error/email_SC_shifted.m: compares the approximation of subgraph centrality obtained with the shifted Arnoldi approach [1, Sec. 5.2.2] to the 'exact' quantity computed with expm for the temporal E-Mail EU network.
 - Example_6_3_Numerical_approximation_error/email_SCres_bipartite.m: compares the upper bound on resolvent-based subgraph centrality computed by Gauss--Radau quadrature applied to the symmetric bipartite network representation to the 'exact' quantity computed with the backslash for the temporal E-Mail EU network.
 - Example_6_3_Numerical_approximation_error/email_SCres_shifted.m: compares the approximation of resolvent-based subgraph centrality obtained with the shifted Arnoldi approach [1, Sec. 5.2.2] to the 'exact' quantity computed with the backslash for the temporal E-Mail EU network.
 - Example_6_3_Numerical_approximation_error/email_TC.m: compares the approximation of total communicability obtained with Krylov subspace methods to the 'exact' quantity computed with expm for the temporal E-Mail EU network.
 - Example_6_3_Numerical_approximation_error/scotland_yard_KC.m: compares the approximation of Katz centrality obtained with Krylov subspace methods to the 'exact' quantity computed with the backslash for the weighted Scotland Yard network.
 - Example_6_3_Numerical_approximation_error/scotland_yard_SC.m: compares Gauss quadrature bounds on subgraph centrality to the 'exact' quantity computed with expm for the weighted Scotland Yard network.
 - Example_6_3_Numerical_approximation_error/scotland_yard_SCres.m: compares Gauss quadrature bounds on resolvent-based subgraph centrality to the 'exact' quantity computed with the backslash for the weighted Scotland Yard network.
 - Example_6_3_Numerical_approximation_error/scotland_yard_TC.m: compares the approximation of total communicability obtained with Krylov subspace methods to the 'exact' quantity computed with expm for the weighted Scotland Yard network.
 - Example_6_4_European_airlines_network/Estrada_index_Gauss_quadrature.m: computes lower and upper bounds on the Estrada index of the European airlines multiplex network [2,3] by applying Gauss, Gauss--Radau, and Gauss--Lobatto quadrature rules to all diagonal entries.
 - Example_6_4_European_airlines_network/KC_MNC_varying_omega.m: computes approximations to marginal layer and marginal node Katz centralities of the European airlines multiplex network [2,3] with the coupling parameter \omega varying between 10^(-2) and 10^3.
 - Example_6_4_European_airlines_network/top_joint_centralities_TC_SC_KC_SCres.m: computes approximations to the top 10 joint centralities of total communicability (TC), subgraph centrality (SC), Katz centrality (KC), and resolvend-based subgraph centrality (SCres) of the European airlines multiplex network [2,3] using Krylov subspace methods and upper Gauss--Radau bounds, respectively.
 - Example_6_4_European_airlines_network/trace_diagonal_estimation_Hadamard.m: computes deterministic estimates of the trace and the diagonal of f(A)=exp(beta*A) using Hadamard vetors in order to obtain approximations to the Estrada index and joint subgraph centralities of the European airlines multiplex network [2,3].
 - Example_6_4_European_airlines_network/trace_diagonal_estimation_Rademacher.m: computes stochastic estimates of the trace and the diagonal of f(A)=exp(beta*A) using stochastic Rademacher vetors in order to obtain approximations to the Estrada index and joint subgraph centralities of the European airlines multiplex network [2,3]. The results for different values of s (number of Rademacher vectors) are averaged over several draws of Rademacher vectors.
 - Example_6_5_Temporal_IMDb_Sci-Fi_network/diagonal_estimation_SC_Hadamard.m: computes deterministic estimates of the diagonal of f(A)=exp(beta*A) using Hadamard vetors in order to obtain approximations to joint subgraph centralities of the temporal IMDb Sci-Fi network.
 - Example_6_5_Temporal_IMDb_Sci-Fi_network/Gauss_quadrature_selected_diagonal_entries.m: computes upper Gauss--Radau quadrature bounds for a specified number of node-layer pairs top-ranked in diagonal estimation.
 - Example_6_5_Temporal_IMDb_Sci-Fi_network/marginal_centralities_TC_KC.m: computes approximations of broadcaster and receiver total communicability and Katz centrality via Krylov subspace methods for the temporal IMDb Sci-Fi network.

**Functions:**
 - subroutines/gauss_lobatto_resolvent.m: computes upper Gauss--Lobatto quadrature bound on resolvent-based subgraph centrality.
 - subroutines/gauss_lobatto_subgraph.m: computes upper Gauss--Lobatto quadrature bound on subgraph centrality.
 - subroutines/gauss_radau_resolvent.m: computes lower or upper Gauss--Radau quadrature bound on resolvent-based subgraph centrality.
 - subroutines/gauss_radau_subgraph.m: computes lower or upper Gauss--Radau quadrature bound on subgraph centrality.
 - subroutines/gauss_resolvent.m: computes lower Gauss quadrature bound on resolvent-based subgraph centrality.
 - subroutines/gauss_subgraph.m: computes lower Gauss quadrature bound on subgraph centrality.
 - subroutines/lanczos_tridiag.m: implementation of the symmetric Lanczos method [7].
 - subroutines/lanczos_tridiag_fh.m: implementation of the symmetric Lanczos method [7] if matrix-vector products are available via function handle.
 - subroutines/MV_prod.m: function that performs the matrix-vector product A*b where A is the supra-adjacency matrix A = blkdiag(A_blkdiag{1}, ... , A_blkdiag{L}) + omega * kron(A_tilde, ones(n,n)) when only A_blkdiag and A_tilde are given, i.e., without explicitly forming A.
 - subroutines/MV_prod_bipartite_adjacency.m: function handle that performs the matrix-vector product [0, A; A', 0]*b with the symmetric bipartite network representation where A is the supra-adjacency matrix.
 - subroutines/MV_prod_transposed.m: function that performs the matrix-vector product A'*b where A is the supra-adjacency matrix.
 - subroutines/funm_kryl/*: funm_kryl toolbox for Matlab [4], which realizes restarted Lanczos and Arnoldi methods for the approximation of f(A)b, see the description in subroutines/funm_kryl/index.htm.
 
 **Version:**
 
 All codes were tested on Matlab R2020b on Ubuntu 20.04.3 LTS as well as on Matlab R2021b on Ubuntu 18.04.6 LTS.
 
 **References:**
 
 - [1] Bergermann, K. & Stoll, M. (2021) Fast computation of matrix function-based centrality measures for layer-coupled multiplex networks. https://arxiv.org/abs/2104.14368v3
 - [2] Cardillo, A., Gómez-Gardenes, J., Zanin, M., Romance, M., Papo, D., Del Pozo, F. & Boccaletti, S. (2013) Emergence of network features from multiplexity. Sci. Rep., 3(1), 1-6. https://doi.org/10.1038/srep01344
 - [3] Taylor, D. (2021) Code Release: Supracentrality. Available at https://github.com/taylordr/Supracentrality
 - [4] Güttel, S. (2008) funm_kryl toolbox for Matlab. Available at http://guettel.com/funm_kryl
 - [5] Paranjape, A., Benson, A. R. & Leskovec, J. (2017) Motifs in temporal networks. In Proceedings of the Tenth ACM International Conference on Web Search and Data Mining, pages 601–610. https://doi.org/10.1145/3018661.3018731
 - [6] Fenu, C. & Higham, D.J. (2017) Block matrix formulations for evolving networks, SIAM Journal on Matrix Analysis and Applications, vol. 38, no. 2, pp. 343-360. https://doi.org/10.1137/16M1076988
 - [7] Golub, G. H. & Van Loan, C. F. (2013) Matrix Computations, volume 3. JHU press, USA.
 
**Authors:**
 - Kai Bergermann, TU Chemnitz, kai.bergermann(at)math.tu-chemnitz.de
 - Martin Stoll, TU Chemnitz, martin.stoll(at)math.tu-chemnitz.de
 
