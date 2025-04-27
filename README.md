# lsgeo: Latent Space Geometry Analysis in Network Models

**lsgeo** is an R package designed to help researchers explore and characterize the latent geometries underlying network models.  
It offers tools for estimating curvature, conducting bootstrap-based hypothesis testing, and classifying networks into Euclidean, spherical, or hyperbolic spaces.

---

## Installation

To get started, you can install the package directly from GitHub:

```r
# Install devtools if needed
install.packages("devtools")

# Install lsgeo
devtools::install_github("jitongj/LS_Geo", subdir = "lsgeo")
```

---

## Overview

Understanding the latent structure of a network can provide valuable insights into its formation mechanisms and global properties.  
In many models, the existence of a latent Euclidean, spherical, or hyperbolic space influences the likelihood of connections between nodes.

The **lsgeo** package implements procedures for:
- Testing competing hypotheses about latent geometries using a bootstrap-based approach.
- Estimating curvature parameters when non-Euclidean geometries are suspected.
- Simulating networks embedded in different latent spaces.
- Classifying observed networks according to their most likely underlying geometry.

These methods are broadly applicable to biological, social, and technological networks.

---

## Main Functions

| Function | Purpose |
|:---|:---|
| `Return_Classification()` | Classify the network's latent geometry. |
| `Bootstrap_Romano()` | Test the appropriateness of a specified latent geometry. |
| `Estimate_Curvature_Sphere()`, `Estimate_Curvature_Hyperbolic()` | Estimate curvature parameters for spherical or hyperbolic spaces. |
| `Make_Euclidean_D()`, `Make_Separated_D_Spherical_3D()`, `Generate_Hyperbolic()` | Generate synthetic node positions and distance matrices. |
| `GenerateY()`, `findCliques_GivenY()` | Generate adjacency matrices and detect cliques. |
| `ComputePerformance()` | Evaluate the performance of geometry classification through simulations. |

All functions are fully documented with examples to help users quickly adapt them to their own networks.

---

## Repository Structure

- `/lsgeo/`: Core R package code and documentation.
- `/examples/`: Supplementary scripts demonstrating applications, including real and simulated data.

---

## Example Applications

- **Neural Networks**: Classifying the latent space of the C. elegans connectome. (See `examples/Classification_CElegans_Repeat.R`)
- **Simulation Studies**: Assessing type I error control and test power across different network settings. (See `examples/Type1_Power.R`)
