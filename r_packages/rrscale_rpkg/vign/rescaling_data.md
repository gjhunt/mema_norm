In this vignette we'll look at how to use rrscale to re-scale data and help discover latent effects.

First we are going to generate data that is the concatenation of two log-normal groups. We do this by taking the outer product between two i.i.d. log-normal vectors, creating group 1 as:

``` r
set.seed(919)
u1 = rlnorm(100)
v1 = rlnorm(100)
Y1 = u1%*%t(v1)
```

and similarly group 2 as

``` r
u2 = rlnorm(100)
v2 = rlnorm(100)
Y2 = .5+u2%*%t(v2)
```

and then we concatenate these together to make a full data matrix (adding some noise)

``` r
Y_nn = rbind(Y1,Y2)
Y = Y_nn + array(rlnorm(prod(dim(Y_nn)),0,.05),dim(Y_nn))
```

Notice that its difficult to tell the groups apart:

``` r
library('reshape2')
library('ggplot2')
group = factor(rep(c(1,2),each=100))
levels(group) = c("group1","group2")
mY = melt(data.frame(Y,group),id.vars="group")
ggplot(data=mY,mapping=aes(x=value,color=group))+geom_histogram(bins=100)+geom_vline(data=aggregate(value~group,data=mY,mean),mapping=aes(xintercept=value,linetype=group),size=1.5)
```

<img src="rescaling_data_files/figure-markdown_github/unnamed-chunk-4-1.png" width="100%" />

Indeed if we look a t-test between the row means across groups we see no difference

``` r
t.test(rowMeans(Y)[group=="group1"],rowMeans(Y)[group=="group2"])
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  rowMeans(Y)[group == "group1"] and rowMeans(Y)[group == "group2"]
    ## t = -0.57743, df = 197.89, p-value = 0.5643
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -1.0335057  0.5653475
    ## sample estimates:
    ## mean of x mean of y 
    ##  3.896319  4.130398

Let's try this after transforming the data with rrscale:

``` r
library('rrscale')
scl = rrscale(Y,run_parallel=FALSE)
```

after running this we get an estimated transformation to help recover the latent group effect. The element "T\_name" tells us that the best transformation is a box-cox-like transformation

``` r
scl$T_name
```

    ## [1] "box_cox_negative"

and the element "par\_hat" tells us the optimal value for the parameter to this transformation:

``` r
scl$par_hat
```

    ## [1] -0.6987309

we can grab the pre-computed RR transformation from the call to rrscale

``` r
trans_Y = scl$RR
```

or we can use the returned "rr\_fn" to calcluate this transformation, they are identical

``` r
trans_Y2 = scl$rr_fn(Y)
all(trans_Y2==trans_Y,na.rm=TRUE)
```

    ## [1] TRUE

Notice that if we plot the transformed Y we see that the group difference is easier to see:

``` r
tmY = melt(data.frame(trans_Y,group),id.vars="group")
ggplot(data=tmY,mapping=aes(x=value,color=group))+geom_histogram(bins=100)+geom_vline(data=aggregate(value~group,data=tmY,mean),mapping=aes(xintercept=value,linetype=group),size=1.5)
```

<img src="rescaling_data_files/figure-markdown_github/unnamed-chunk-11-1.png" width="100%" /> indeed the t-test is now significant

``` r
t.test(rowMeans(trans_Y)[group=="group1"],rowMeans(trans_Y)[group=="group2"])
```

    ## 
    ##  Welch Two Sample t-test
    ## 
    ## data:  rowMeans(trans_Y)[group == "group1"] and rowMeans(trans_Y)[group == "group2"]
    ## t = -3.297, df = 182.51, p-value = 0.001175
    ## alternative hypothesis: true difference in means is not equal to 0
    ## 95 percent confidence interval:
    ##  -0.5156955 -0.1295520
    ## sample estimates:
    ##  mean of x  mean of y 
    ## -0.1613588  0.1612650

If we plot he first two PCs for the transformed and un-transformed data we can see the group difference much better after transformation:

``` r
plot(svdc(Y)$u[,1:2],col=group)
```

<img src="rescaling_data_files/figure-markdown_github/unnamed-chunk-13-1.png" width="100%" />

``` r
plot(svdc(trans_Y)$u[,1:2],col=group)
```

<img src="rescaling_data_files/figure-markdown_github/unnamed-chunk-14-1.png" width="100%" />

Here we are using the "svdc" function from the rrscale package which calculates 'completed" right and left singular vectors in the presence of missing values. We can also look at the canonical correlation between the group and the first two PCS for the transformed and untransformed data:

``` r
cancor(model.matrix(~1+group),svdc(Y)$u[,1:2])
```

    ## $cor
    ## [1] 0.6986683
    ## 
    ## $xcoef
    ##                   [,1]
    ## groupgroup2 -0.1414214
    ## 
    ## $ycoef
    ##             [,1]        [,2]
    ## [1,] -0.07388689 -1.58366207
    ## [2,]  0.99745825 -0.09317112
    ## 
    ## $xcenter
    ## (Intercept) groupgroup2 
    ##         1.0         0.5 
    ## 
    ## $ycenter
    ## [1] -0.054830536 -0.002675599

``` r
cancor(model.matrix(~1+group),svdc(trans_Y)$u[,1:2])
```

    ## $cor
    ## [1] 0.9767767
    ## 
    ## $xcoef
    ##                  [,1]
    ## groupgroup2 0.1414214
    ## 
    ## $ycoef
    ##            [,1]      [,2]
    ## [1,] -0.4454289 0.8979825
    ## [2,]  1.6282824 0.6951998
    ## 
    ## $xcenter
    ## (Intercept) groupgroup2 
    ##         1.0         0.5 
    ## 
    ## $ycenter
    ## [1]  0.002759017 -0.058307115

and we can see that it is much higher for the transformed data signifiying these principal components capture the latent group structure better.
