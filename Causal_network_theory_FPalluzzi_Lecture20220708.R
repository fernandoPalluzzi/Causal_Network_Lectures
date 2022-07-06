


        #-----------------------------------------------#
        #                                               #
        #  PART I - Structural Equation Modeling (SEM)  #
        #                                               #
        #-----------------------------------------------#



# 1.1. Start R and change to your current working directory.



# 1.2. Data inspection and SEM fitting.


# SPECT brain data for controls (group = 1), Frontotemporal dementia 
# (FTD) patients carrying Granulin gene mutations (GRN+; group = 2), and 
# FTD Granulin-wt (GRN-) patients (group = 3).

# SPECT (Single Photon Emission Computed Tomography) is a brain imaging 
# technique used to detect altered blood flow.

# Dataset from: Premi E, et al. 2013 (PMID: 23687363)

# Brain areas chart
#    Temp: temporal
#  Pariet: parietal
#   DLPFC: dorsolateral prefrontal cortex
#   Orbit: orbitofrontal
# AntCing: anterior cingulate

brain.data <- read.table("FTD_GRN_brain_SPECT.csv",
                         header = TRUE,
                         sep = ";",
                         dec = ",")

head(brain.data)


# Load required libraries

library(SEMgraph)
library(SEMdata)


# Build the directed interaction network as a data.frame

#         from        to    weight
# 1     Temp.d   Pariet.d        1
# 2     Temp.d    DLPFC.d        1
# 3     Temp.d  AntCing.d        1
# 4   Pariet.d    Orbit.d        1
# 5    Orbit.d  AntCing.d        1
# 6    DLPFC.d    Orbit.d        1
# 7    DLPFC.d  AntCing.d        1

G <- data.frame(from = c("Temp.d", "Temp.d", "Temp.d", "Pariet.d", "Orbit.d", "DLPFC.d", "DLPFC.d"),
                to = c("Pariet.d", "DLPFC.d", "AntCing.d", "Orbit.d", "AntCing.d", "Orbit.d", "AntCing.d"),
                weight = rep(1, 7))


# Convert the dataframe into an igraph network object

G <- graph_from_data_frame(G, directed = TRUE)

plot(G)    # igraph

gplot(G)   # SEMgraph


# Some aesthetics ...

V(G)
V(G)$color <- c("yellow", "lightblue", "yellow", "yellow", "yellow")

E(G)
E(G)$color <- c("black", "red3", "black", "black", "red3", "red3", "black")
E(G)$width <- c(1, 2, 1, 1, 2, 2, 1)

G

gplot(G, fontsize = 24)


# Create the path diagram (causal network) corresponding to the brain 
# interaction network

path.diagram <- '
Pariet.d ~ Temp.d
DLPFC.d ~ Temp.d
Orbit.d ~ DLPFC.d + Pariet.d
AntCing.d ~ Temp.d + DLPFC.d + Orbit.d'


# Extract SPECT data only (no ID, no group)

names(brain.data)

semdata <- brain.data[, 3:14]

head(semdata)


# Data inspection

# Data standardization (centering and scaling).
# Note: 
# SEM effect estimates are based on correlation and standardization 
# does not affect it.

semdata <- apply(semdata, 2, scale)

dim(semdata)


# Correlation matrix

semcor <- cor(semdata)
semcor

# Correlation matrix determinant
# The determinant should be > 0. A very small (i.e., close to zero) 
# value could be a symptom of multicollinearity.
# Multicollinearity happens when one predictor in a multiple regression 
# can be linearly predicted from other predictors, with high accuracy.

det(semcor)


# Let us reduce the dataset to the model-implied variables

predictors <- V(G)$name

semdata <- semdata[, colnames(semdata) %in% predictors]

semcor <- cor(semdata)
det(semcor)


# SEM fitting (common model)

fit <- sem(path.diagram, data = brain.data, fixed.x = TRUE, se = "robust.huber.white")

summary(fit, rsquare = TRUE, fit.measures = TRUE)



# 1.3. SEM fitting interpretation.


# 1.3.a. Parameter estimates.

# Effect sizes, covariances, and variances are model parameter estimates.

# The effect size (often referred to as "beta"), or direct effect (DE), 
# can be seen as the magnitude of the directed causal effect x -> y, 
# from one source x to its target node y.
# In terms of linear regression, the magnitude of the DE can be interpreted 
# as the variation induced in y for a unit increase in x, by keeping 
# constant all the other regressors acting on y.
# In R, this relationship is indicated by: y ~ x
# (i.e., the endogenous varible io always on the left side of the ~).

# The covariance between two variables can be viewed as a bidirected 
# interaction x <-> y between two variables, and represents the joint 
# variability between them. The Pearson's (linear) correlation coefficient 
# is the normalized version of the covariance, corresponding to the 
# ratio between the cov(x, y) and the product of their standard deviations.
# In R (SEM notation) the covariance is represented by: x ~~ y.

# Finally, the variance is a measue of dispersion of the values of a 
# random variable from their average value. If a variable y is explained 
# in terms of a predictor f(X1, X2, ..., Xm), the proportion of variance 
# of y that is not explained by f(.) is referred to as residual variance.
# Therefore, the variance can be intended as a general measure of 
# variability.

# When we deal with a system composed by different variables, its 
# variance-covariance matrix represents the co-variation structure of 
# the system itself. 

# When these model parameters are estimated, given data and an estimation 
# criterion (e.g., maximum likelihood estimation), we also give an error 
# measure for that estimation (e.g., the standard error).
# In a linear regression model, the standardized effect is represented 
# by its score z = (b - beta0)/SE, where SE is the standard error and 
# beta0 (usually, equal to 0) is the effect size under null hypothesis.
# If data distribution deviates from normality (this is often the case), 
# robust SE estimation methods, including Huber-White "sandwitch" or 
# bootstrap estimation, can be use to reduce the impact of normality 
# constraints violation.

# The SE allow us to calculate also the 95% confidence interval (CI) 
# and P-value for each effect size estimate:

# CI95 = (beta - 1.96*SE, beta + 1.96*SE);
# where 1.96 is the z value at a significance level of 0.05;

# pvalue = 2*pnorm(-abs(z)); for a two-tailed test.

# The 95% CI is the range of estimates for a given (unknown) parameter.
# If this interval contains the beta0 value (i.e., the beta value under 
# null hypothesis), then the direct effect is non-significant.
# Note that this definition of significance does not require P-value 
# calculation.


# 1.3.b. Fitting performance (chi-squared test and fitting indices).

# In general, a model "well fits data" when it is capable to explain 
# large part of data variability, usually expressed in terms of 
# covariance structure. In any case, the better the model goodness-of-fit, 
# the higher the likelihood that the observed data has been generated from 
# that model (or an equivalent one).


# Chi-squared (Model Test User Model)

# This chi-squared test evaluates the divergence of the current model 
# from a saturated model (i.e., completely connected network), under
# the null hypothesis H0: Sigma(theta) = E(S);
# i.e., there is no difference between the model covariance matrix 
# Sigma(theta) and the expected value of the sample covariance matrix E(S).

# A non-significant test will not reject H0, and the current model 
# will be not significantly different from a saturated one: our model 
# explains (almost) all the variability explained by the completely 
# connected one.


# RMSEA (Root Mean Squared Error of Approximation)

# Unfortunately, the chi-squared statistic is biased towards high 
# dimensionality: it tends to reject the null hypothesis for very small 
# differences between the current model and the saturated one.
# The RMSEA provides a sample size-adjusted version of the chi-squared 
# statistic. An RMSEA < 0.05 and non-significant RMSEA P-value denotes 
# a good model fitting.
# Reference:
# Byrne, B.M. Structural Equation Modeling With AMOS: Basic Concepts, 
# Applications, and Programming, Second Edition. Routledge (2013).


# SRMR (Standardized Root Mean squared-Residual)

# This fitting index directly computes the discrepancy between the model 
# covariance matrix and its expected value E(S). The smaller the SRMR 
# the better the fitting.
# Generally, an SRMR < 0.08 denotes a good model fitting.
# Reference:
# Hu, L. & Bentler, P.M. Cutoff Criteria for Fit Indexes in Covariance 
# Structure Analysis: Conventional Criteria Versus New Alternatives. 
# Structural Equation Modeling, 6(1), 1-55 (1999).


# CFI (Comparative Fit Index)

# Differently from the RMSEA and SRMR indices, the CFI compares the 
# current model with a baseline one. Sometimes, RMSEA and SRMR are 
# defined as "badness-of-fit" indices, while the CFI is the actual GOF 
# index. A CFI >= 0.95 denotes a good model fitting.
# Reference:
# Bentler, P.M. Comparative fit indexes in structural models. 
# Psychol Bull 107, 238–246 (1990).


# 1.3.c. R-squared.

# Within the system of regression equations, the outcome of each equation 
# is an endogenous variable. In general, every variable receiving at least 
# an incoming connection is an endogenous variable.
# Variables with only outgoing connections are instead exogenous variables.
# The R-squared is the proportion of variance of an endogenous variable 
# that is explained by its predictors (i.e., its incoming connections).
# The higher the R-squared values of the endogenous variables the higher 
# the proprotion of variance explained by the whole model.



# 1.4. Evaluating system perturbation and indirect effects (IE).


# Outcome recoding: FTD (1), healthy (0)

table(brain.data$groups)

y <- ifelse(brain.data$groups == 1, 0, 1)


# Using SEMgraph::SEMrun function to evaluate the diesease-induced perturbation

sem <- SEMrun(G, semdata, y)

# Group perturbation effect over system variables:
# the variation induced in a variable by the phenoptype, 
# moving from the control group to the cases one.

# Evaluating system status
summary(sem$fit, rsquare = TRUE, fit.measures = TRUE)

gplot(sem$graph, fontsize = 24)


# Adding more parameters (IE and covariance)

path.ie <- '

# Path diagram
DLPFC.d ~ b1*Temp.d
Pariet.d ~ b2*Temp.d
Orbit.d ~ b3*Pariet.d + b4*DLPFC.d
AntCing.d ~ b5*Temp.d + b6*Orbit.d + b7*DLPFC.d

# Covariance
Pariet.d ~~ DLPFC.d

# Indirect effect
IE1:= b1*b4*b6
IE2:= b1*b7'

#  ~ direct effect (DE)
# ~~ lavaan syntax to define a new covariance
# := lavaan syntax to define a new model parameter other than DE or covariance

# Model fitting

fit <- sem(path.ie, data = brain.data, fixed.x = TRUE, se = "robust.huber.white")

summary(fit, rsquare = TRUE, fit.measures = TRUE)



#----------------------------------------------------------------------#



        #-----------------------------------------------#
        #                                               #
        #   PART II - Molecular pathways and networks   #
        #                                               #
        #-----------------------------------------------#



# Start R and change to your current working directory.

# Loading required libraries:

library(SEMgraph)
library(SEMdata)



# 2.1. Exploring an igraph object: interaction networks.

kegg
reactome
string

length(V(kegg))

length(E(kegg))

V(kegg)$name

E(kegg)$weight

table(E(kegg)$weight)



# 2.2. Molecular pathways.

kegg.pathways

list(kegg.pathways)

length(kegg.pathways)

head(names(kegg.pathways), 30)
head(names(reactome.pathways), 30)

pathway.list <- names(kegg.pathways)

j <- grep("Amyotrophic", pathway.list, fixed = TRUE)
j
pathway.list[j]

als <- kegg.pathways[[j]]
als

gplot(als, l = "circo")



# 2.3. Exploring biological networks.


# Graph connected components

als.components <- properties(als)

G1 <- als.components[[1]]

G2 <- als.components[[2]]


# Creating a "label" attribute for gene symbols.
# Gene symbols could be ambiguous and contain special characters.
# For this reason, we use entrez gene IDs as node names.
# If labels are defined, they will be visualized in the graph plot 
# instead of node names.

library(org.Hs.eg.db)

V(G2)$label <- mapIds(org.Hs.eg.db, V(G2)$name, column = "SYMBOL", keytype = "ENTREZID")

gplot(G2)
gplot(G2, l = "circo")



# 2.4. Testing normality for expression (count) data.


# Let's make an example with the C9orf72 gene.

# Entrez ID for C9orf72 gene: "203228".
# See https://www.genecards.org/ for details.

x <- alsData$exprs[, colnames(alsData$exprs) == "203228"]

# Log2 transform is often used to reduce the "heavy-tail" effect of
# biological count data, due to overdispersion.

# WARNING: THIS DATASET IS ALREADY LOG2-TRANSFORMED!
# The following line of code is for didactic purposes only.
# x <- log2(x + 1)
# The "+ 1" is referred to as pseudocount. This is needed to prevent the 
# logarithm argument from going to 0.

# Let's simulate the Gaussian distribution with same mean and SD 
# as the C9orf72 gene expression distribution.

x.gauss <- rnorm(n = 10000000, mean = mean(x), sd = sd(x))

# Plot of the reference Gaussian and the C9orf72 expression distributions

plot(density(x.gauss), lwd = 3, col = "darkblue", lty = 2)
lines(density(x), lwd = 3, col = "red3")

# Assessing normality: the quantile-quantile (q-q) plot

library(car)
qqPlot(x, distribution = "norm")

# Assessing normality: the Shapiro-Wilk test

H <- shapiro.test(x)
H

# Applying the Shapiro-Wilk test in batch

p <- apply(alsData$exprs, 2, function(x) shapiro.test(x)$p.value)

# How many genes show non-Gaussian behaviour?

length(p)
length(p[p <= 0.05])


# Nonparanormal transorm

library(huge)

als.npn <- huge.npn(alsData$exprs)

p.npn <- apply(als.npn, 2, function(x) shapiro.test(x)$p.value)

length(p.npn[p.npn <= 0.05])



# 2.4. Testing system perturbation in ALS.


# Let's see how ALS perturbs the G2 component of the network.

G2.sem <- SEMrun(G2, als.npn, alsData$group)

G2 <- G2.sem$graph
V(G2)$label <- mapIds(org.Hs.eg.db, V(G2)$name, column = "SYMBOL", keytype = "ENTREZID")

gplot(G2, l = "circo")


# Extracting differentially regulated nodes (DRNs) and edges (DREs)

G1.sem <- SEMrun(G1, als.npn, alsData$group)

DRN <- G1.sem$gest[G1.sem$gest$pvalue < 0.05,]
head(DRN, 10)

G1.sem2 <- SEMrun(G1, als.npn, alsData$group, fit = 2)

DRE <- G1.sem2$dest[G1.sem2$dest$pvalue < 0.05,]
head(DRE, 10)



# 2.5. Estimation of a new model from data.


# Let's create a unique large directed interaction network from the 
# fusion of KEGG and Reactome ones

reference <- graph.union(kegg, reactome)
reference


# Four different model estimation strategies (use ?modelSearch for deails):


#  BASIC. No reference network needed. The input network is only used to 
#         define the initial topological order.

# model <- modelSearch(G1, als.npn, search = "basic", beta = 0)


# DIRECT. The reference network is used to validate new inferred connections.
#         Only new direct connections are allowed.

# model <- modelSearch(G1, als.npn, gnet = reference, search = "direct", beta = 0)


#  INNER. The reference network is used to validate new inferred connections.
#         Mediated connections are allowed, but new directed paths must 
#         involve mediators from the input network.
#         A maximum of d - 1 mediators is allowed.

# model <- modelSearch(G1, als.npn, gnet = reference, d = 2, search = "inner", beta = 0)


# OUTER. The reference network is used to validate new inferred connections.
#        Mediated connections are allowed.
#        Mediators can be imported from the reference network.
#        A maximum of d - 1 mediators is allowed.

model <- modelSearch(G1, als.npn, gnet = reference, d = 2, search = "outer", beta = 0)

# The beta argument defines the threshold for the minimum absolute LASSO 
# beta coefficient. If beta = 0 (default), all the non-zero estimated 
# connections will be included within the input network.

# Let's call G our new graph model
G <- model$graph

V(G)$label <- mapIds(org.Hs.eg.db, V(G)$name, column = "SYMBOL", keytype = "ENTREZID")
gplot(G)


# 2.6. Defining network seeds.


# Using DRNs as seeds. Seeds are nodes with special features (e.g., they 
# are differentially regulated or known for being associated with the 
# phenotype of interst).
# Often, one is interested in weighting them (e.g.: 1 for "seed" and 
# 0 for "non-seed") to give them importance when applying search/extraction 
# algorithm, including tree search, random walks, or clustering.

sem <- SEMrun(G, als.npn, alsData$group)

DRN <- sem$gest[sem$gest$pvalue < 0.05,]
head(DRN)

# Seed list
DRN <- rownames(DRN)
DRN



# 2.7. Edge weighting and automatied seed detection.


# SEMgraph offers three ways for weighting network edges
# (use ?weightGraph for details):

# 1. "r2z"
#    Fisher's r-to-z transform (r = Pearson's correlation).
#    Weights are sign and P-value for the corresponding z-score.
#    This is the default method (faster and more interpretable).

# 2. "sem"
#    Weights are sign and P-value of z = w/SE(w), where w is a new 
#    parameter defined as the weighted sum of the total effect of the 
#    group on source and sink nodes, adjusted by node degree centrality.

# 3. "cov"
#    Weights are sign and P-value of z = w/SE(w), where w is a new 
#    parameter combining the group effect on the source node (mean group 
#    difference, adjusted by source degree centrality), the sink node 
#    (mean group difference, adjusted by sink degree centrality), and 
#    the source-sink interaction (correlation difference). 

G <- weightGraph(G, als.npn, alsData$group)

summary(E(G)$pv)

table(E(G)$zsign)

# If one does not have a criterion for seed detection, but still wants 
# to define them, SEMgraph yields three sets of seeds, if the "seed" 
# argument is enabled:

G <- weightGraph(G, als.npn, alsData$group, seed = c(0.05, 0.5, 0.9))

# The seed argument include three cutoffs. Firstly, the significance 
# level of the direct group effect over each node. Secondly, the 
# prototype clustering distance (= 1 - |r|) cutoff. Finally, 
# the closeness percentile (0.5 is the second quartile, aka the median).

table(V(G)$pvlm)

table(V(G)$proto)

table(V(G)$qi)



# 2.8. Applying network reduction algorithms over weighted networks.


# Three-based solution: the Kou implementation of the Steiner tree problem.

R <- activeModule(G, type = "kou", seed = DRN, eweight = "pvalue")

V(R)$label <- mapIds(org.Hs.eg.db, V(R)$name, column = "SYMBOL", keytype = "ENTREZID")
gplot(R)

R <- SEMrun(R, als.npn, alsData$group)$graph

V(R)$label <- mapIds(org.Hs.eg.db, V(R)$name, column = "SYMBOL", keytype = "ENTREZID")
gplot(R)

r <- properties(R)

# The Steiner tree is a minimum-cost/maximum-information flow graph.
# It maximizes the information flow while minimizing the overall cost 
# (often defined as the sum of network distances), without cycles.
# If the initial graph is directed (with no bi-directed edges) the 
# resulting tree is a directed acyclic graph (DAG).
# This representation is useful in molecular biology to represent 
# signaling cascades (e.g., from the membrane receptor to the nucleous,
# to modulate transcriptional regulation).


# Random walkers: when the network has a modular structure.

W <- SEMdci(G, als.npn, alsData$group, type = "wtc", method = "BH", alpha = 0.05)

W <- SEMrun(W, als.npn, alsData$group)$graph

V(W)$label <- mapIds(org.Hs.eg.db, V(W)$name, column = "SYMBOL", keytype = "ENTREZID")
gplot(W)

w <- properties(W)

gplot(w[[1]])
gplot(w[[1]], l = "circo")

# A random walk (RW) is the stochastic process of taking random steps 
# within a mathematical space: in our case, the graph W(V, E).
# If W has a modular structure is likely that the RW will get stuck 
# inside a module (i.e., community or cluster) after a number of steps.
# In many RW-based algorithms a potential dissipation function establish 
# the length of a certain RW, allowing us to detect network modules.


# Average causal effect (ACE): finding multiple cascades and directed graphlets.

A <- SEMdci(G, als.npn, alsData$group, type = "ace", method = "BH", alpha = 0.05)

A <- SEMrun(A, als.npn, alsData$group)$graph

V(A)$label <- mapIds(org.Hs.eg.db, V(A)$name, column = "SYMBOL", keytype = "ENTREZID")
gplot(A)

a <- properties(A)

# In its general definition, the total effect (TE) between two nodes 
# x and y in a path P = x -> ... -> y is the sum of the DE x -> y 
# (i.e., the "beta" coefficient) and the IE = b(1,2)*b(2,3)*...*b(k-1,k).
# If there is no direct effect, TE = IE, while in absednce of indirect 
# effect(s) TE = DE. Therefore, the TE is useful to assess the presence 
# of a causal effect between any pair of connected nodes in a directed 
# network.
# One convenient way of estimating the TE is through the definition of 
# ACE (Pearl J, 1998; https://doi.org/10.1177/0049124198027002004).
# The simplest estimation of the TE as ACE is possible in DAGs, thorugh 
# linear regression. The parent set pa(X) of X blocks all backdoor 
# (i.e., confounding) paths from X to Y, and the ACE is equal to the
# b(Y,X|Z) coefficient in a multiple regression Y ~ X + pa(X).



#----------------------------------------------------------------------#
