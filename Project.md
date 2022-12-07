# CX 4640 Final Project
## Understanding Singular Value Decomposition(SVD) and its Applications
---
###### Name: Niv Magril
###### Topic: 11
###### Title: The SVD. What is it, what it can be used for, how it can be computed.
---
# Table of Contents
1. [Quick History](#history)
2. [Base Idea and Introduction](#introduction)
3. [Constructing SVD and Example](#example)
4. [Numerical Analysis Applications of SVD](#analysis)
    - [Matrix Properties](#matrices)
    - [Pseudoinverse](#pseudoinverse)
    - [Lower-Rank Matrix Approximation](#lower)
    - [Total Least Squares](#total)
6. [Real World Applications](#world)
    - [Image Processing and Compression](#image)
    - [Facial Recognition](#face)
    - [Quantitative Finance](#quant)
    - [Web Searching](#web)
8. [Coding SVD](#coding)
9. [References](#references)


### Quick History<a name="history"></a>


For those of us who enjoy an origin story here is a quick timeline of events regarding the development of SVD. We start with Eugenio Beltrami and Camille Jordan who found the SVD for simplification of bilinear forms in 1870s. Jordan went on and obtained geometric interpretation of the largest singular value. Next J. J. Sylvester wrote two papers on the SVD in 1889. In these papers he found algorithms to diagonalise quadratic and bilinear forms by means of orthogonal substitutions. This was followed by Erhard Schmidt who in 1907 discovered the SVD for function spaces while investigating integral equations. The main problem he was investigating was finding the best rank k approximations to $A$ of the form, $u_{1}v_{1}^{t} + \cdots + u_{k}v_{k}^{t}$<sup>[[10]](#Hestenes)</sup>.

After a few years we had Autonne who found the SVD for complex matrices in 1913. This was succeeded by Eckhart and Young who managed to extend SVD to rectangular matrices in 1936<sup>[[3]](#young)</sup>. Then we witnessed Golub and Kahan introducing SVD in numerical analysis in 1965, the main reason we highlight SVD in this presentation<sup>[[7]](#Kahan)</sup>. Finally, Golub expanded on his numerical analysis insights by proposing an algorithm for SVD in 1970<sup>[[8]](#Golub)</sup>.


### Base Idea and Introduction<a name="introduction"></a>
Singular Value Decomposition is as simple as its name, its a method of decomposing matrix into matrices composed of singular values. Of course it is a little more intuitvie than that. Given a matrix, $A \in \mathbb{R}^{m \times n}$, we can factorize this matrix into the following form,
$$A = U\Sigma V^{T}$$
These matrices are defined below:
- U: An $m \times m$ unitary matrix 
  - Values $u_{i}$ are known as the left singular vectors \& $UU^{\ast} = I$
- V: An $n \times n$ unitary matrix 
  -  Values $v_{i}$ are known as the left singular vectors \& $VV^{\ast} = I$
- $\Sigma$: A $m \times n$ diagonal matrix with non-negative real numbers on the diagonal
  - Where the entries along the diagonal are denoted $\sigma_{i}$ and are usually order in decreasing order: $\sigma_{1} \geq \sigma_{2} \geq \sigma_{3} \geq \dots$

$\lceil$ As a quick aside we define a matrix to be unitary if it can be matrix multiplied with its respective hermetian matrix and results in the identity matrix. That is for a matrix $X \in \mathbb{R}^{m \times n}$ we have that $XX^{H} = X^{H}X = I.$ More information can be found at the following [link](https://en.wikipedia.org/wiki/Unitary_matrix), but we can continue with our introduction of SVD. $\rfloor$

As you can see, SVD decomposes the matrix into 3 different matrices. Of these 3 matrices we have that $U$ and $V$ are unitary matrices such that the following holds,
$$UU^{T} = U^{T}U = I$$
$$VV^{T} = V^{T}V = I$$

These properties prevents us from working with a possible "ugly matrix" by producing matrices that are easier to handle as illustrated above. However, what is SVD actually doing? So you have a matrix A, which is the matrix you want to decompose using SVD. This is a transformation matrix that transforms a group of vectors to new space. Let’s define the original orthogonal vectors as $v$'s and the transformed orthogonal vectors as $u$'s. At last, let’s normalize the transformed matrix so that it becomes easier to handle the results. As a matter of fact, this is what SVD is doing. It’s basically dividing different transformations into each matrix $U$, $\Sigma$, and $V$. A more geometrical explanation of SVD can be found [here](https://gregorygundersen.com/blog/2018/12/10/svd/). We will now go into how we produce these 3 matrices.









### Constructing SVD and Example<a name="example"></a>

We note that in general there is **no exact method for computing the SVD**. However, we will now illustrate the one basic algorithm that produce the matrices listed above and use a guided example to see how it works. 

First we will produce the singular values $\sigma_{i}$ of our diagonal matrix $\Sigma$. Given a matrix $A$ the singula values are computed from the eigenvalues of $AA^{T}$. We see an example of this below,

$$A = \begin{bmatrix}
3 & 2 & 2\\
2 & 3 & -2
\end{bmatrix}\ \ \ AA^{T} = \begin{bmatrix}
17 & 8\\
8 & 17
\end{bmatrix}$$

Now we have that the eignevalues of this matrix are computed from the following,
$det(AA^{T} - \lambda I) = \lambda^{2} - 34\lambda + 225 = (\lambda - 25)(\lambda - 9)$.

From this we have that our eigenvalues and subsequent singular values are given as follows $\sigma_{1} = \sqrt{25} = 5$ and $\sigma_{2} = \sqrt{9} = 3$, and thus we have the following.

$$ \Sigma  = \begin{bmatrix}
5 & 0 & 0\\
0 & 3 & 0
\end{bmatrix} $$

Now we find the right singular vectors, the columns of V, by finding an orthonormal set of eigenvectors of $A^{T}A$. Now the eigenvalues of $A^{T}A$ are also given by the eigenvalues of $AA^{T}$. Now since $AA^{T}$ is symmetrix we know that the eigenvectors will be orthogonal. So now we compute the eigenvector for $\lambda_{1} = 25$,

$$A^{T}A - \lambda_{1}I = \begin{bmatrix}
-12 & 12 & 2\\
12 & -12 & -2\\
2 & -2 & -17\\
\end{bmatrix}$$ 

Reducing this matrix will give us our first eigenvector which is given by,

$$v_{1} = \begin{bmatrix}
1/\sqrt{2}\\
1/\sqrt{2}\\
0\end{bmatrix}$$

Similarly we perform the same computation for $\lambda_{2} = 9$.

$$A^{T}A - \lambda_{2}I = \begin{bmatrix}
4 & 12 & 2\\
12 & 4 & -2\\
2 & -2 & -1\end{bmatrix}$$ 

Reducing this matrix will give us our first eigenvector which is given by,

$$v_{2} = \begin{bmatrix}
1/\sqrt{18}\\
-1/\sqrt{18}\\
4/\sqrt{18}
\end{bmatrix}$$

Now for the last eigenvector we need a unit vector that is perpendicular or orthogonal to both $v_{1}$ and $v_{2}$. In this case to be perpendicular to 

$$v_{1} = \begin{bmatrix} 
a\\
b\\
c \end{bmatrix}$$

we need $a = -b$. Moreover, we need $v_{2}^{T}v_{3} = 0$. From these conditions we have that,

$$v_{3} = \begin{bmatrix}
a\\
-a\\
-a/2
\end{bmatrix} = \begin{bmatrix}
2/3\\
-2/3\\
-1/3
\end{bmatrix}$$

Thus we have all of our orthogonal eigenvectors and we can define the following,

$$V = \begin{bmatrix}
1/\sqrt{2} & 1/\sqrt{18} & 2/3\\
1/\sqrt{2} & -1/\sqrt{18} & -2/3\\
0 & 4/\sqrt{18} & -1/3\\
\end{bmatrix}$$

Finally we can compute $U$ through the formula $\sigma u_{i} = Av_{i}$ or $u_{i} = \frac{1}{\sigma}Av_{i}$. Doing so we can compute $U$ and we can find our full Singular Value Decomposition below,

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





### Numerical Analysis Applications of SVD<a name="analysis"></a>

#### Matrix Properties<a name="matrices"></a>
Below we note a few important applications of SVD in regards to linear algebra and numberical analysis. One main benefit is its ability to compute other aspects regarding matrices that we list below<sup>[[9]](#text)</sup>:
- **Minimum norm solution:** given a system of equations $Ax\cong b$ we can determine the following,
$$x = \sum_{\sigma_{i}\neq0} \frac{u_{i}^{T}b}{\sigma_{i}}v_{i}$$
- **Euclidean matrix norm:** Can be given by $||A|| _{2} = \sigma_{max}$
- **Euclidean condition number of matrix:** $cond_{2}(A) = \frac{\sigma_{max}}{\sigma_{min}}$
- **Rank of a matrix:** rank(A) = number of nonzero singular value, count( $\sigma_{i} \neq 0$)
- **Range of a matrix:** Given by the left singular vectors of $U$ corresponding to non-zero singular values
- **Null space of a matrix:** Given by the right singular vectors of $V$ corresponding to the zeroed singular values.

#### Pseudoinverse<a name="pseudoinverse"></a>

Additionally the SVD can be used to find the **pseudoinverse** of our matrix $A$, this is given by 
$$A^{+} = V\Sigma^{+}U^{T}$$

Where $\Sigma^{+}$ is the pseudoinverse of our diagonal matrix, that is each pseudoinverse of our scalar $\sigma_{i}$ is defined to be $\frac{1}{\sigma_{i}}$ for $\sigma_{i}\neq0$.

We note that the $A^{+}$ exists regardless whether or not our matrix is square or has full rank. Additionaly, if $A$ is square and nonsingular we have that $A^{+}=A^{-1}.

The Pseudoinverse is very useful in computing the minimum norm solution, that is $Ax\cong b$ is given by $x=A^{+}b

#### Lower-Rank Matrix Approximation<a name="lower"></a>

Suppose we want to best approximate a matrix A by a rank-k matrix. This is where the SVD comes in handy. Our singular value decomposition can also be written in the form


$$A = U\Sigma V^{T} = \sigma_{1}E_{1} + \sigma_{2}E_{2} + \cdots + \sigma_{n}E_{n},$$

where $E_{i} = u_{i}v_{i}^{T}$. Now we see that each $E_{i}$ has a rank of 1 and can be stored using only $m+n$ storage locations. Thus the product $E_{i}x$ can be computed using only $m+n$ multiplication. Condensed approximations of A are obtained by omitting from summation terms corresponding to small singular values.


Below we show a step by step process of computing this rank-k approximation:
1. Compute the SVD of our matrix $A = U\Sigma V^{T}$
2. Keep only the top $k$ right singular vectors: set $V_{k}^{T}$ equal to the first $k$ rows of our matrix $V^{T}$
3. Keep only the top $k$ left singular vectors: set $U_{k}$ equal to the first $k$ columns of U
4. Keep only the top $k$ singular values: set $\Sigma_{k}$ equal to the first $k$ rows and columns of $\Sigma$, corresponding to the $k$ largest singular values of $A$.
5. The rank-k approximation is then given by $A_{k} = U_{k}\Sigma_{k}V_{k}^{T}$


These approximation come in handy for image processing, data compression, information retreival, and more topics we will go over in the following section.

#### Total Least Squares<a name="total"></a>

In general least squares is applicable when our right hand side, b, is subject to some random error, but our matrix $A$ is known accurately. However, in the case that $A$ is also subject to error we can apply total least square. 

Computing the SVD of $[A, b]$ gives us the total least squares solution by minimizing the orthogonal distance between model and data. The solution turns out to be the right-singular vector, $v_{m}$, of $A$ corresponding to the smallest singular value.

More specificaly, singular value decomposition can be used to find a unique solution to total least squares problems.



### Real World Applications<a name="world"></a>


Singular Value Decomposition has been not been hoarded only for the domain of numerical analysis but is applicable in multiple other fields. We illustrate some of these real world applications below

#### Image Processing and Compression<a name="image"></a>

 Here we consider how SVD can be used to produce reduced image sizes. We begin by understanding that large images are formed by correspondingly large matrices, requiring alot of memory to store the image. By rewriting the image in its pixel by pixel form and removing the smaller singular values, we can form smaller matrices which would require less memory storage. We would lose some refinement with each loss of a singular value, but we would retain the overall image feature. An example of this can be seen below.
<div align="center">

<img src="https://github.com/magriln/MATH-4640-Project/blob/82fdedb40f049934835c19eb3b51f8d99f39ff8e/GrayscaleSVD.png"  width="400" height="200">
  
  
<b>Figure 1</b>: Illustrating Image Compression by removal of Singular Values and Subsequent Pixels
</div>

  
SVD can also be applied to the saturation of each pixel. Of course this becomes more nuanced when evaluating colored images, but we can also see how the saturation of each pixel can also affect grayscale images. An example of the effect of SVD on saturation can be seen below.
  
<div align="center">
  
<img src="https://github.com/magriln/MATH-4640-Project/blob/82fdedb40f049934835c19eb3b51f8d99f39ff8e/ImageSaturation.png"  width="400" height="200">

  
<b>Figure 2</b>: Analyzing the Effect of SVD on Grayscale Saturation</div>
  
Finally within the domain of grayscale images we can see how the number of pixels and their respective saturation can be combined to illustrate the full effect of SVD on image compression.
  
<div align="center">
  
  
<img src="https://github.com/magriln/MATH-4640-Project/blob/82fdedb40f049934835c19eb3b51f8d99f39ff8e/ImageValues.png"  width="600" height="200">
 
<b>Figure 3</b>: Comparing Singular Value Removal in Regards to Pixel Count and Saturation
</div>

  
The process of SVD can also be expanded to full color images. Each pixel in full color image has color saturation representation values of 0 to 255 for RGB. This entails complexity to the image, which requires a greater amount of memory to store an image. Representing each color relative to the full color image, we are able to see the amount of contribution each color has. In order to implement the SVD process we have to first separate the full color image into its red, green, and blue layers, as each of these three colors has its own matrix of information for the image. We remove the smallest singular values from each of the color matrices, and then we reconstruct the full color image using the modified color matrices.

<div align="center">
  
<img src="https://github.com/magriln/MATH-4640-Project/blob/82fdedb40f049934835c19eb3b51f8d99f39ff8e/ColoredSVD.png"  width="400" height="400">

<b>Figure 4</b>: SVD Being Applied to Color Images</div>

  
By applying the process of Singular Value Decomposition to images by using pixel saturation matrices for grayscale or full color images, we can compress the storage size while retaining the number of pixels. Isolating the least important pieces of information that are stored in the images and have removed them methodically, leaving only the most important components of the images. This process of removing the smallest singular values from the saturation matrices allows us to retain as much of the image quality as possible.

#### Facial Recognition<a name="face"></a>

Over the past few decades numerous facial recognition algorithms have been explored. Under different lighting conditions, poses and facial expressions we have witnessed progress towards recognition. A facial recognition algorithm and its implementation can be considered as a system similar to the image compression we highlighted above. Inputting a two dimensional images, a predetermined library of faces distinguishes the input image as a user’s face. With the hopeful output being the ability to discern face images.

Since the facial recognition problem itself deals with images, we treat this similar to image compression. That is each pixel in an image is considered as a coordinate in an n-dimensional space, where n is the total number of pixels per image. Now given a collection of images we form a library of faces that serve as contenders for recognition. Since each image is a point in n-space, it would be computationally efficient to reduce the overall storage space of each image. 

Using the eigenface technique, we form the space of images, which are then projected into a low dimensional space using singular value decomposition. Thus the high dimensional n-space is transformed to a set of uncorrelated singular values that span most if not all variation in the original data set. The determined singular values can be used to thereby attempt to reconstruct an input image, and subsequently
classify a face as a library element. We can see an example of the construction of eigenface below.

<div align="center">
  
<img src="https://github.com/magriln/MATH-4640-Project/blob/e76c9e723fa967c2991ce74aee9591c2e9d4fb84/eigenfaces.jpg"  width="600" height="200">

<b>Figure 5</b>: Illustrating the Deconstruction of Facial Images Using SVD</div>

From these images we can then input similar faces and compress them to identify their similarities to our library of faces. COmparing singular values we can correctly identify the faces of complete strangers or at least identify their defining features.


#### Quantitative Finance<a name="quant"></a>

We see SVD being used in domains including finances but it plays a core piece in financial modeling. For example stock prices are affected by multiple factors, and various methods have been proposed to improve prediction accuracy. However, not all of the proposed features are valid, and there is often noise in the features—such as political, economic, and legal factors—which can lead to poor prediction results. Using SVD we can reconstruct the features of stock data, eliminate data noise, retain the most effective data features, and improve the accuracy of prediction.  We can reconstruct the data by selecting the large singular values, which will reduce if not eleiminate the noise in the data and improve data quality. 

In stocks we can construct a matrix matrix that is $m$ by $n$, where m is the number of stock data records and n is the number of stock features. Then applying neural network algorithms such as LSTM are applied to predict the behavior of these stocks. An example of this relationship is illustrated in the image below.

<div align="center">
  
<img src="https://github.com/magriln/MATH-4640-Project/blob/827eeefa620f18bfa287ce00ebc16a2d9645a588/StocksSVD.png"  width="400" height="300">

<b>Figure 6</b>: Illustrating the Use of SVD-LSTM in Stock Predictions</div>

#### Web Searching<a name="web"></a>

Search engines like Google use enormous matrices of cross referencing checking what words are on each page. Upon a Google search, the higher ranks of this matrix usually go to pages with your key words that have lots of links to them. But there are billions of pages out there, and storing a billion by billion matrix is trouble. This is not considering querying through it.

SVD comes in handy in this application as well. In searching, we really only care about the main directions that the Web is taking, the top results. So the first few singular values create a very good approximation for the enormous matrix, can be searched relatively quickly and provide compression ratios of millions to one. 



Of course these are just a handfull of real world applications, and we can see Singular Value Decomposition being applied everywhere from medical fields to athletics. We will now analyze how most of these industries code their own SVD by illustrating a base algorithm.


### Coding SVD<a name="coding"></a>

As mentioned in the history portion of this paper Golub and Kahan constructed a universal algorithm for constructing the SVD. This was further altered for special cases of SVD. We illustrate the pseudocode for Golub-Kahan and Golub-Reinsch<sup>[[2]](#algorithms)</sup> algorithm below,

<div align="center">
  
<img src="https://github.com/magriln/MATH-4640-Project/blob/e3215f104a70e65e6656ccf66a99ad670a9f1e40/SVD1.png"  width="450" height="300">
<img src="https://github.com/magriln/MATH-4640-Project/blob/e3215f104a70e65e6656ccf66a99ad670a9f1e40/SVD2.png"  width="350" height="300">

<b>Figure 7</b>: Illustrating Golub-Kahan, and Golub-Reinsch SVD Algorithms</div>



---
### References<a name="references"></a>
---

1. Algebraic structures. Linear Algebra Homepage. (n.d.). Retrieved December 7, 2022, from http://staff.imsa.edu/~fogel/LinAlg/ 
2. Cline, A. K., &amp; Dhillon, I. S. (2021). Computation of the decomposition - university of Texas at Austin. https://www.cs.utexas.edu/. Retrieved December 7, 2022, from https://www.cs.utexas.edu/~inderjit/public_papers/HLA_SVD.pdf <a name="algorithms"></a>
3. Dey, S. (2018, January 8). Eigenfaces and a simple face detector with PCA/SVD in python. sandipanweb. Retrieved December 7, 2022, from https://sandipanweb.wordpress.com/2018/01/06/eigenfaces-and-a-simple-face-detector-with-pca-svd-in-python/ 
4. Eckart, C.; Young, G. (1936). "The approximation of one matrix by another of lower rank". Psychometrika. 1 (3): 211–8. doi:10.1007/BF02288367. S2CID 10163399.<a name="young"></a>
5. El Ghaoui, L. (2021). Applications of SVD. Retrieved December 7, 2022, from https://inst.eecs.berkeley.edu/~ee127/sp21/livebook/l_svd_apps.html 
6. Ernstberger, S. L. (2020). Singular Value Decomposition: Applications to Image Processing . Undergraduate Research. Retrieved December 7, 2022, from https://www.lagrange.edu/academics/undergraduate/undergraduate-research/index.html 
7. Gavin, H. P. (2017, December 17). Total least squares - duke university. Retrieved December 7, 2022, from https://people.duke.edu/~hpgavin/SystemID/CourseNotes/TotalLeastSquares.pdf 
8. Golub, Gene H.; Kahan, William (1965). "Calculating the singular values and pseudo-inverse of a matrix". Journal of the Society for Industrial and Applied Mathematics, Series B: Numerical Analysis. 2 (2): 205–224. Bibcode:1965SJNA....2..205G. doi:10.1137/0702016. JSTOR 2949777.<a name="Kahan"></a>
9. Golub, G. H.; Reinsch, C. (1970). "Singular value decomposition and least squares solutions". Numerische Mathematik. 14 (5): 403–420. doi:10.1007/BF02163027. MR 1553974. S2CID 123532178.<a name="Golub"></a>
10. Heath, M. T. (2018). Scientific computing: An introductory survey. Society for Industrial and Applied Mathematics. <a name="text"></a>
11. Hestenes, M. R. (1958). "Inversion of Matrices by Biorthogonalization and Related Results". Journal of the Society for Industrial and Applied Mathematics. 6 (1): 51–90. doi:10.1137/0106005. JSTOR 2098862. MR 0092215.<a name="Hestenes"></a>
12. Roughgarden, T., &amp; Valiant, G. (2022, April 24). CS 168: The Modern Algorithmic Toolbox, spring 2022. The Modern Algorithmic Toolbox (CS168), Spring 2022. Retrieved December 7, 2022, from https://web.stanford.edu/class/cs168/index.html 
13. Rutgers. (n.d.). SVD and Signal Processing. Applied Optimum Signal Processing. Retrieved December 7, 2022, from http://eceweb1.rutgers.edu/~orfanidi/aosp/ 
14. Verma, J. K. (2020, March 13). Introduction to linear algebra, fifth edition (2016). Introduction to Linear Algebra, 5th Edition. Retrieved December 7, 2022, from https://math.mit.edu/~gs/linearalgebra/ 
15. Wikimedia Foundation. (2022, November 7). Singular value decomposition. Wikipedia. Retrieved December 7, 2022, from https://en.wikipedia.org/wiki/Singular_value_decomposition#Total_least_squares_minimization 
