===================================================================
README: Datasets

There are eight datasets included in this folder:

1.	Adopted.mat
2.	fiberpaper.dat
3.	irisf.mat
4.	Rohwer.mat
5.	sales.txt
6.	VocabGrowth.mat
7.	waterstrider.mat
8.	wheatprotein.txt

===================================================================

1.	Adopted Children

Description:

     Data are a subset from an observational, longitudinal, study on
     adopted children.  Is child's intelligence related to intelligence
     of the biological mother and the intelligence of the adoptive
     mother?

     The child's intelligence was measured at age 2, 4, 8, and 13 for
     this sample.  How does intelligence change over time, and how are
     these changes related to intelligence of the birth and adoptive
     mother?

Usage:

     load Adopted.mat
     
Format:

     A data frame with 62 observations on the following 6 variables.

Value:

     Column 1: adoptive mother's years of education (proxy for her IQ)

     Column 2: biological mother's score on IQ test

     Column 3: IQ of child at age 2

     Column 4: IQ of child at age 4

     Column 5: IQ of child at age 8

     Column 6: IQ of child at age 13

Source:

     Ramsey, F.L. and Schafer, D.W. (2002). _The Statistical Sleuth: A
     Course in Methods of Data Analysis (2nd ed)_, Duxbury.

     This data set is identical to ‘ex1605’ in the ‘Sleuth2’ package.

References:

     Friendly, Michael (2010). HE Plots for Repeated Measures Designs.
     _Journal of Statistical Software_, 37(4), 1-40. URL <URL:
     http://www.jstatsoft.org/v37/i04/>.

     Skodak, M. and Skeels, H.M. (1949). A Final Follow-up Study of One
     Hundred Adopted Children, _Journal of Genetic Psychology_ *75*:
     85-125.

See Also:

     ‘ex1605’

===================================================================

2.	Pulp and paper property

Description:

     This data set contains measurements of properties of pulp fibers
     and the paper made from them.

Usage:

     load fiberpaper.dat
     
Format:

     A matrix with 62 observations on 8 variables.

Value:

     Column 1: Breaking length.

     Column 2: Elastic modulus.

     Column 3: Stress at failure.

     Column 4: Burst strength.

     Column 5: Arithmetic fiber length.

     Column 6: Long fiber fraction.

     Column 7: Fine fiber fraction.

     Column 8: Zero span tensile.

References:

     Johnson, R.A. and Wichern, D.W. (2007). Applied Multivariate
     Statistical Analysis, 6th edition.

===================================================================

3.	Fisher's iris data

Description:

     This data set contains measurements of 4 characteristics for 3
     species of iris.

Usage:

     load irisf
     
Format:

     A matrix with 150 observation on 5 variables.

Value:

     Column 1: Indicator of the species: 1 for Iris Setaso, 2 for Iris
          Versicolor and 3 for Iris Verginica.

     Column 2: Sepal length.

     Column 3: Sepal width.

     Column 4: Petal length.

     Column 5: Petal width.

References:

     Johnson, R.A. and Wichern, D.W. (2007). Applied Multivariate
     Statistical Analysis, 6th edition.

===================================================================

4. Rohwer Data Set

Description:

     Data from an experiment by William D. Rohwer on kindergarten
     children designed to examine how well performance on a set of
     paired-associate (PA) tasks can predict performance on some
     measures of aptitude and achievement.

Usage:

     load Rohwer.mat

Format:

     A data frame with 69 observations on the following 5 variables.

Value:

     Column 1: ‘SAT’ a numeric vector: score on a Student Achievement Test

     Column 2: ‘PPVT’ a numeric vector: score on the Peabody Picture Vocabulary
          Test

     Column 3: ‘Raven’ a numeric vector: score on the Raven Progressive Matrices
          Test

     Column 4: ‘n’ a numeric vector: performance on a 'named' PA task

     Column 5: ‘s’ a numeric vector: performance on a 'still' PA task

Details:

     The variables ‘SAT’, ‘PPVT’ and ‘Raven’ are responses to be
     potentially explained by performance on the paired-associate (PA)
     learning task‘n’ and ‘s’.

Source:

     Timm, N.H. 1975).  Multivariate Analysis with Applications in
     Education and Psychology.  Wadsworth (Brooks/Cole), Examples 4.3
     (p. 281), 4.7 (p. 313), 4.13 (p. 344).

References:

     Friendly, M. (2007).  HE plots for Multivariate General Linear
     Models.  Journal of Computational and Graphical Statistics,
     *16*(2) 421-444.  <URL: http://datavis.ca/papers/jcgs-heplots.pdf>

===================================================================

5.	The quality of its sales staff

Description:

     A firm is attempting to evaluate the quality of its sales staff
     and is trying to find an examination or series of tests that may
     reveal the potential for good performance in sales.  The firm has
     selected a random sample of 50 sales people and has evaluated each
     on 3 measures of performance: growth of sales, profitability of
     sales, and new-account sales. These measures have been converted
     to a scale, on which 100 indicates "average" performance.  Each of
     the 50 individuals took each of 4 tests, which purported to
     measure creativity, mechanical reasoning, abstract reasoning, and
     mathematical ability, respectively.

Usage:

     load sales.dat
     
Format:

     A matrix with 50 observations on 7 variables.

Value:

     Column 1: Index of sales growth.

     Column 2: Index of sales profitability.

     Column 3: Index of new-account sales.

     Column 4: Score on creativity test.

     Column 5: Score on mechanical reasoning test.

     Column 6: Score on abstract reasoning test.

     Column 7: Score on mathematics test.

References:

     Johnson, R.A. and Wichern, D.W. (2007). Applied Multivariate
     Statistical Analysis, 6th edition.

===================================================================

6.	Vocabulary growth data

Description:

     Data is collected from the Laboratory School of the University of
     Chicago.  The data consists of scores from a cohort of pupils in
     grades 8-11 on the vocabulary section of the Cooperative Reading
     Test.  The scores are scaled to a common, but arbitrary origin and
     unit of measurement, so as to be comparable over the four grades.

Usage:

     load VocabGrowth.mat
     
Format:

     A matrix with 64 observations on 4 variables.

Value:

     Column 1: Grade 8 vocabulary score.

     Column 2: Grade 9 vocabulary score.

     Column 3: Grade 10 vocabulary score.

     Column 4: Grade 11 vocabulary score.

References:

     Bock, R.D. (1975). Multivariate Statistical Methods in Behavioral
     Research.

===================================================================

7.	Water strider data

Description:

     Water striders live on the surface of lakes and ponds.  These
     insects grow in five distinct stages before reaching adulthood,
     these stages are called instars; at each transition thy shed their
     skins, which are also their skeletons.  This data set contains
     measurements of 8 characteristics for 3 species of water striders:
     L. dissortis, L. rufoscutellatus and L. esakii.  The measurements
     are taken for the first three instars, when the sex of the water
     striders are undetermined.

Usage:

     load waterstrider
     
Format:

     A matrix with 90 observations on 10 variables.

Value:

     Column 1: Binary indicator: 1 for L. esakii and 0 for other species.

     Column 2: Binary indicator: 1 for L. dissortis and 0 for other species.

     Column 3: Logarithm (natural base) of length of the first antennal
          segment.  Raw measurements are in millimeters.

     Column 4: Logarithm (natural base) of length of the second antennal
          segment.  Raw measurements are in millimeters.

     Column 5: Logarithm (natural base) of length of the third antennal
          segment.  Raw measurements are in millimeters.

     Column 6: Logarithm (natural base) of length of the fourth antennal
          segment.  Raw measurements are in millimeters.

     Column 7: Logarithm (natural base) of length of middle femora.  Raw
          measurements are in millimeters.

     Column 8: Logarithm (natural base) of length of middle tibiae.  Raw
          measurements are in millimeters.

     Column 9: Logarithm (natural base) of length of hind femora.  Raw
          measurements are in millimeters.

     Column 10: Logarithm (natural base) of length of hind tibiae.  Raw
          measurements are in millimeters.

References:

     Klingenberg, C.R. and Spence, J.R. (1993). Heterochrony and
     Allometry: Lessons from the Water Strider Genus Limnoporus.
     Evolution 47, 1834 - 1853.

===================================================================

8.	The protein content of ground wheat samples.

Description:

     The data are the result of an experiment to calibrate a near
     infrared reflectance (NIR) instrument for measuring the protein
     content of ground wheat samples. The protein content of each
     sample (in percent) was measured by the standard Kjeldahl method.
     In Fearn (1983), the problem is to find a linear combination of
     the measurements that predicts protein content. The estimated
     coefficients can then be entered into the instrument allowing the
     protein content of future samples to be read directly. The first
     24 cases were used for calibration and the last 26 samples were
     used for prediction.

Usage:

     load wheatprotein.dat
     
Format:

     A matrix with 50 observations on 8 variables.

Value:

     Column 1 to Column 6: Measurements of the reflectance of NIR radiation
          by the wheat samples at 6 wavelength in the range 1680-2310
          nm.  The measurements were made on the log(1/reflectance)
          scale.

     Column 7: The protein content of each sample (in percent).

     Column 8: Binary indicator, 0 for high protein content and 1 for low
          protein content. The cut off point is if the protein content
          is smaller than 9.75.

References:

     Fearn, T. (1983). A misuse of ridge regression in the calibration
     of a near infrared reflectance instrument.

