## (R)obust (R)e-scaling

The Robust re-scaling transformation (RR) is a transformation the help reveal latent structure in data. It uses three steps to transform the data:

1. Gaussianize the data,
2. z-score Transform the data,
3. remove extreme outliers from the data.

The sequence of these transformations helps focus analyses on consequential variance in the data rather than  having it be focused on variation resulting from the feature's  measurement scale or outliers.

This is the github page for our R package 'rrscale' to perform the transformation.

A basic vignette using rrscale is available [here](vign/rescaling_data.md) ([html](vign/rescaling_data.html) [rmd](vign/rescaling_data.Rmd))

