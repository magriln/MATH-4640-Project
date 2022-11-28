# CX 4640 Final Project
## Understanding Singular Value Decomposition(SVD) and its Applications
---
###### Name: Niv Magril
###### Topic: 11
###### Title: The SVD. What is it, what it can be used for, how it can be computed.
----
### Base Idea and Introduction
###### Singular Value Decomposition is as simple as its name, its a method of decomposing matrix into matrices composed of singular values. Of course it is a little more intuitvie than that. Given a matrix, $A \in \mathbb{R}^{m \times n}$ we can factorize this matrix into the following form,
$$A = U\Sigma V^{T}$$
###### These matrices are defined below:
- U: An $m \times m$ orthogonal matrix 
  - Values $u_{i}$ are known as the left singular vectors 
- V: An $n \times n$ orthogonal matrix 
  -  Values $v_{i}$ are known as the left singular vectors
- $\Sigma$: A $m \times n$ diagonal matrix with non-negative real numbers on the diagonal
  - Where the entries along the diagonal are denoted $\sigma_{i}$ and are usually order in decreasing order: $\sigma_{1} \geq \sigma_{2} \geq \sigma_{3} \geq \dots$

###### 
