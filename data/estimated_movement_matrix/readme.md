## Question 
Assuming a fully connected network with directed and weighted air transport links (flows), each node represents one country, I have the total volumn of inflow and outflow of each node in this network, how to estimate the directed flows for each link between two nodes. Please note that I only have the total inflow and outflow for every node, but I don't know the individual flow from / to each of the other nodes, which is the quantity I want to estimate. Please also note that each node does not have flow to itself. Give the solution in R.

## Gravity model (for ICAO data)
To estimate the directed flows for each link between two nodes in a fully connected network with directed and weighted air transport links, you can use the Gravity Model. The Gravity Model is commonly used to predict the flow between nodes (i.e., countries) based on certain attributes like total flow capacity, distance, and economic factors.

Here's a basic outline of how to apply the Gravity Model in this scenario:

1. **Gather necessary data:**
    - Total flow for every node
    - Distance between each pair of nodes
    - Any other relevant variables or information that affect the flow (e.g., GDP, population, etc.)

2. **Define the Gravity Model equation:**

    The general equation for the Gravity Model is:

    ```
    T_ij = k * (O_i ^ alpha) * (D_j ^ beta) / D_ij ^ gamma
    ```
    
    Where:
    - `T_ij` represents the estimated flow from node i to node j.
    - `O_i` represents the outflow volume of node i.
    - `D_j` represents the inflow volume of node j.
    - `D_ij` represents the distance between node i and node j.
    - `alpha`, `beta`, and `gamma` are scaling parameters to adjust the weights of outflows, inflows, and distances, respectively.
    - `k` is a proportionality constant.

3. **Calibrate the model parameters:**

    To calibrate the model parameters `alpha`, `beta`, `gamma`, and `k`, you'll need additional data on actual flow volumes between some pairs of nodes. If you don't have any, you may start with some initial values (e.g., `alpha = 1`, `beta = 1`, `gamma = 2`) and see how the model performs, then adjust the parameters accordingly.

4. **Calculate the estimated flow between nodes:**

    With the calibrated parameters, use the Gravity Model equation to estimate the directed flows for each link between two nodes. Make sure you normalize the flow values so that they satisfy the constraint of total inflow and outflow at every node.

**Details of the gravity model in R** please refer to https://www.ajs.or.at/index.php/ajs/article/view/688

### GLM and gravity model
The generalized linear model (GLM) is a flexible generalization of ordinary linear regression that allows for response variables with error distribution models other than a normal distribution. In the case of a gravity model, we use GLM to capture the relationships between multiple predictor variables and the dependent variable while accommodating the non-linear and non-normal nature of the data.

A gravity model usually aims to describe the interaction or flow between spatial units, such as countries, states, or cities. The main idea behind such a model is that the interaction between two locations is proportional to their "mass" (e.g., population, GDP) and inversely proportional to the distance between them.

In our example, we used GLM with a Poisson distribution (commonly used for count data) and a log link function. By incorporating a log link function, we can specify multiplicative relationships among the predictors. This aligns well with the gravity model's basic assumptions:

The origin's mass (GDP per capita and population) has a positive effect on the interaction.
The destination's mass (GDP per capita and population) also has a positive effect on the interaction.
The distance between the two locations has a negative effect on the interaction.
With this configuration, the GLM estimates coefficients for each predictor in the gravity model, considering their log-transformed form. Once the coefficients are estimated, you can make predictions for traffic flows or other interaction measures between spatial units. The combination of GLM and gravity model concepts allows for more accurate analysis and prediction of interactions or flows between spatial units.

## Pairwise interaction (fir WTO data)
In the context of this problem, pairwise interaction represents the relative strength or intensity of interactions between each pair of countries (nodes) based on their normalized inflows and outflows. The interaction value gives an initial estimate for how much flow is expected to exist between two countries.

To calculate pairwise interaction, we use an outer product operation between the normalized outflows and normalized inflows. It is computed such that:

`pairwise_interaction[i, j] = normalized_outflows[i] * normalized_inflows[j]`

This generates a square matrix where the number of rows and columns equals the number of countries (nodes). Each cell in the matrix `[i, j]` gives the preliminary interaction between country `i` and country `j`. In essence, it assumes that the intensity of interaction between two countries depends on both their normalized inflow and outflow values.

However, these initial pairwise interaction values do not guarantee that the sum of flows leaving a particular node will match its original outflow value, which is required in the final solution. To adjust for this, we renormalize the interaction matrix row-wise such that the row sums match the original outflow values, while preserving the relative strengths of pairwise interactions among all nodes.

### Relationship to basic gravity model
The pairwise interaction model described in the previous response is, in fact, a variation of the basic gravity model. The gravity model in its most fundamental form is used to predict the interactions between two entities based on their size and distance. In the context of international trade or migration, it predicts the amount of flow between two countries based on their economic sizes (such as GDP) and geographical distances.

The basic gravity model is given by:

`flow[i, j] = G * size[i] * size[j] / distance[i, j]^p`

Where:
- `G` is a constant.
- `size[i]` and `size[j]` represent the 'size' of country `i` and country `j` respectively, which can be economic size (e.g., GDP), population, etc.
- `distance[i, j]` represents the geographical distance between the centroids of the two countries.
- `p` is an exponent, usually close to 2, which represents the impact of distance on the flow.

In the pairwise interaction model discussed earlier, we have normalized inflows and outflows instead of size factors and omitted the distance factor. Thus, the pairwise interaction model focuses purely on the interaction strength based on normalized inflows and outflows:

`pairwise_interaction[i, j] = normalized_outflows[i] * normalized_inflows[j]`

As you can see, the pairwise interaction model only considers the inflow and outflow values for calculating the flow between countries. It omits the distance factor and any other size-related variables, making it a simplified gravity model tailored specifically for problems involving flows among nodes where distance or size is not considered.

