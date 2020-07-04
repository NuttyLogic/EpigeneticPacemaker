# Epigenetic Pacemaker Algorithm

## EPM Description

Given $i$ methylation sites and $j$ individuals a single methylation site can be described
as $\hat{m}_{ij} = m^0_i + r_is_j + \epsilon _{ij}$ where $\hat{m}_{ij}$ is the observed methylation value, $m^0_i$ is
the initial methylation values, $r_i$ is the rate of change, $s_j$ is the epigenetic state, and $\epsilon _{ij}$
is a normally distributed error term. Given an input matrix $\hat{M} = [\hat{m_{ij}}]$ the goal of the EPM is
find the optimal values of $m^0_i$, $r_i$, and $s_j$ to minimize the error between the predicted and
observed methylation values across a system of methylation sites. Under the EPM $m^0_i$ and $r_i$ are characteristic
of the site for all individuals and $s_j$ is shared by all sites within a system of methylation sites for every
individual.

The EPM optimization is accomplished through an implementation of a fast conditional expectation maximization algorithm
that maximizes the model likelihood by minimizing the residual sum of squares error. When fitting the EPM each
methylation site is assigned an independent rate of change and starting methylation value, while each individual is
assigned an epigenetic state. The initial epigenetic state is provided by the user and should represent a best guess.
The epigenetic state is then updated through each iteration of the EPM to  minimize the error across
the observed epigenetic landscape. Because the $s_j$ is updated while fitting the EPM the condition of linearity
between the methylation values and trait of interest is relaxed.

## EPM Implementation

The EPM algorithm is implemented as follows

1. fit $i$ site models using the user provided state predictions to get $r_i$ and $m_0$
2. update $s_j$ to minimize $\epsilon_{ij}^2$
      - $s_j = \frac{\sum_{i \leq n} r_i(\hat{m_{ij}} - m^0_i)}{\sum_{i \leq n} r^2_i}$
3. refit site models using $s_j$
4. repeat step 2 and 3 until model improvements $\leq$ specified threshold or maximum number of iterations reached
