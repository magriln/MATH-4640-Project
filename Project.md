# CX 4640 Final Project
## Understanding Singular Value Decomposition(SVD) and its Applications
---
###### Name: Niv Magril
###### Topic: 11
###### Title: The SVD. What is it, what it can be used for, how it can be computed.
----
### Base Idea and Introduction
###### Singular Value Decomposition is as simple as its name, its a method of decomposing matrix into matrices composed of singular values. Of course it is a little more intuitvie than that. Given a matrix, $A \in \mathbb{R}^{m \times n}$, we can factorize this matrix into the following form,
$$A = U\Sigma V^{T}$$
###### These matrices are defined below:
- U: An $m \times m$ unitary matrix 
  - Values $u_{i}$ are known as the left singular vectors 
- V: An $n \times n$ unitary matrix 
  -  Values $v_{i}$ are known as the left singular vectors
- $\Sigma$: A $m \times n$ diagonal matrix with non-negative real numbers on the diagonal
  - Where the entries along the diagonal are denoted $\sigma_{i}$ and are usually order in decreasing order: $\sigma_{1} \geq \sigma_{2} \geq \sigma_{3} \geq \dots$

###### As a quick aside we define a matrix to be unitary when we matrix multiply a matrix and its respective hermetian matrix to produce the identity matrix. That is for a matrix $X \in \mathbb{R}^{m \times n}$ we have that $XX^{H} = X^{H}X = I.$ More information can be found at the following [link](https://en.wikipedia.org/wiki/Unitary_matrix), but we can continue with our introduction of SVD.
###### As you can see, SVD decomposes the matrix into 3 different matrices. Of these 3 matrices we have that $U$ and $V$ are unitary matrices such that the following holds,
$$UU^{T} = U^{T}U = I$$
$$VV^{T} = V^{T}V = I$$
###### These properties prevents us from working with a possible "ugly matrix" by producing matrices that are easier to handle as illustrated above. However, what is SVD actually doing? So you have a matrix A, which is the matrix you want to decompose using SVD. This is a transformation matrix that transforms a group of vectors to new space. Let’s define the original orthogonal vectors as $v$'s and the transformed orthogonal vectors as $u$'s. At last, let’s normalize the transformed matrix so that it becomes easier to handle the results. As a matter of fact, this is what SVD is doing. It’s basically dividing different transformations into each matrix $U$, $\Sigma$, and $V$. We will now go into how we produce these 3 matrices.

### Constructing SVD and Example
