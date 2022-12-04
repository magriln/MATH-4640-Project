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
  - Values $u_{i}$ are known as the left singular vectors \& UU^{\ast} = I
- V: An $n \times n$ unitary matrix 
  -  Values $v_{i}$ are known as the left singular vectors \& VV^{\ast} = I
- $\Sigma$: A $m \times n$ diagonal matrix with non-negative real numbers on the diagonal
  - Where the entries along the diagonal are denoted $\sigma_{i}$ and are usually order in decreasing order: $\sigma_{1} \geq \sigma_{2} \geq \sigma_{3} \geq \dots$

###### As a quick aside we define a matrix to be unitary when we matrix multiply a matrix and its respective hermetian matrix to produce the identity matrix. That is for a matrix $X \in \mathbb{R}^{m \times n}$ we have that $XX^{H} = X^{H}X = I.$ More information can be found at the following [link](https://en.wikipedia.org/wiki/Unitary_matrix), but we can continue with our introduction of SVD.
###### As you can see, SVD decomposes the matrix into 3 different matrices. Of these 3 matrices we have that $U$ and $V$ are unitary matrices such that the following holds,
$$UU^{T} = U^{T}U = I$$
$$VV^{T} = V^{T}V = I$$
###### These properties prevents us from working with a possible "ugly matrix" by producing matrices that are easier to handle as illustrated above. However, what is SVD actually doing? So you have a matrix A, which is the matrix you want to decompose using SVD. This is a transformation matrix that transforms a group of vectors to new space. Let’s define the original orthogonal vectors as $v$'s and the transformed orthogonal vectors as $u$'s. At last, let’s normalize the transformed matrix so that it becomes easier to handle the results. As a matter of fact, this is what SVD is doing. It’s basically dividing different transformations into each matrix $U$, $\Sigma$, and $V$. A more geometrical explanation of SVD can be found [here](https://gregorygundersen.com/blog/2018/12/10/svd/). We will now go into how we produce these 3 matrices.

### Constructing SVD and Example
###### We will now illustrate the basic algorithm to produce the matrices listed above and use a guided example to see how it works. First we will produce the singular values $\sigma_{i}$ of our diagonal matrix $\Sigma$. Given a matrix $A$ the singula values are computed from the eigenvalues of $AA^{T}$. We see an example of this below,
$$A = \begin{bmatrix}
3 & 2 & 2\\
2 & 3 & -2
\end{bmatrix}\ \ \ AA^{T} = \begin{bmatrix}
17 & 8\\
8 & 17
\end{bmatrix}$$

###### Now we have that the eignevalues of this matrix are computed from the following, $det(AA^{T} - \lambda I) = \lambda^{2} - 34\lambda + 225 = (\lambda - 25)(\lambda - 9)$.
###### From this we have that our eigenvalues and subsequent singular values are given as follows $\sigma_{1} = \sqrt{25} = 5$ and $\sigma_{2} = \sqrt{9} = 3$, and thus we have the following.
$$\Sigma  = \begin{bmatrix}
5 & 0 & 0\\
0 & 3 & 0
\end{bmatrix}$$
###### Now we find the right singular vectors, the columns of V, by finding an orthonormal set of eigenvectors of $A^{T}A$. Now the eigenvalues of $A^{T}A$ are also given by the eigenvalues of $AA^{T}$. Now since $AA^{T}$ is symmetrix we know that the eigenvectors will be orthogonal. So now we compute the eigenvector for $\lambda_{1} = 25$,

$$A^{T}A - \lambda_{1}I = \begin{bmatrix}
-12 & 12 & 2\\
12 & -12 & -2\\
2 & -2 & -17\\
\end{bmatrix}$$ 
###### Reducing this matrix will give us our first eigenvector which is given by,
$$v_{1} = \begin{bmatrix}
1/\sqrt{2}\\
1/\sqrt{2}\\
0
\end{bmatrix}$$
###### Similarly we perform the same computation for $\lambda_{2} = 9$.
$$A^{T}A - \lambda_{2}I = \begin{bmatrix}
4 & 12 & 2\\
12 & 4 & -2\\
2 & -2 & -1\\
\end{bmatrix}$$ 
###### Reducing this matrix will give us our first eigenvector which is given by,
$$v_{2} = \begin{bmatrix}
1/\sqrt{18}\\
-1/\sqrt{18}\\
4/\sqrt{18}
\end{bmatrix}$$
###### Now for the last eigenvector we need a unit vector that is perpendicular or orthogonal to both $v_{1}$ and $v_{2}$. In this case to be perpendicular to 
$$v_{1} = \begin{bmatrix} 
a\\
b\\
c 
\end{bmatrix}$$
###### we need $a = -b$. Moreover, we need $v_{2}^{T}v_{3} = 0$. From these conditions we have that,
$$v_{3} = \begin{bmatrix}
a\\
-a\\
-a/2
\end{bmatrix} = \begin{bmatrix}
2/3\\
-2/3\\
-1/3
\end{bmatrix}$$
###### Thus we have all of our orthogonal eigenvectors and we can define the following,
$$
V = \begin{bmatrix}
1/\sqrt{2} & 1/\sqrt{18} & 2/3\\
1/\sqrt{2} & -1/\sqrt{18} & -2/3\\
0 & 4/\sqrt{18} & -1/3\\
\end{bmatrix}$$
###### Finally we can compute $U$ through the formula $\sigma u_{i} = Av_{i}$ or $u_{i} = \frac{1}{\sigma}Av_{i}$. Doing so we can compute $U$ and we can find our full Singular Value Decomposition below,
$$A = U\Sigma V^{T} = \begin{bmatrix}
1/\sqrt{2} & 1/\sqrt{2}\\
1/\sqrt{2} & -1/\sqrt{2}
\end{bmatrix}\begin{bmatrix}
5 & 0 & 0\\
0 & 3 & 0
\end{bmatrix}\begin{bmatrix}
1/\sqrt{2} & 1/\sqrt{2} & 0\\
1/\sqrt{18} & -1/\sqrt{18} & 4/\sqrt{18}\\
2/3 & -2/3 & -1/3\\
\end{bmatrix}$$

### Edge Cases and Faults

###### In general every real matrix has an SVD with real entries, every complex matrix has an SVD with complex entries. However, there can exist SVD implementation that fail to converge









###### Now that we've illustrated the computation and special cases of SVD we will now look into its importance and application.

### Numerical Analysis Applications of SVD
#### Matrix Properties
###### Now from our lectures we noted a few important applications of SVD. One main benefit is its ability to compute other aspects regarding matrices that we list below:
- **Minimum norm solution:** given a system of equations $Ax\cong b$ we can determine the following,
$$x = \sum_{\sigma_{i}\neq0} \frac{u_{i}^{T}b}{\sigma_{i}}v_{i}$$
- **Euclidean matrix norm:** Can be given by $||A|| _{2} = \sigma_{max}$
- **Euclidean condition number of matrix:** $cond_{2}(A) = \frac{\sigma_{max}}{\sigma_{min}}
- **Rank of a matrix:** rank(A) = number of nonzero singular value, count($\sigma_{i}\neq0)

#### Pseudoinverse
###### Additionally the SVD can be used to find the **pseudoinverse** of our matrix $A$, this is given by 
$$A^{+} = V\Sigma^{+}U^{T}$$

