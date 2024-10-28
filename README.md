# Rec

Implement instruction:

For simulation I-IV

For experiment I and II
1. Run exp1data.R to get simulated data sets
2. Run cv2n.R, cv3n.R to get result of CV and choose parameter for splines
3. Run exp1n.R, exp1nn.R, exp1nnn.R to get estimation results of n=400, 600, 800
4. Run exp1nIBS.R, exp1nCindex.R, exp1logrank.R, exp1MAE.R to get prediction results of IBS, C-index, logrank test and MAE
5. Run result.R to analysis estimation results from exp1n.R, exp1nn.R, exp1nnn.R

For experiment III
1. Run exp1data.R to get simulated data sets
2. Run cv2n.R, cv3n.R to get result of CV and choose parameter for splines
3. Run exp1n.R, exp1nn.R, exp1nnn.R to get estimation results of n=400, 600, 800
4. Run test1n.R, test1nn.R, test1nnn.R to get result of Cramer test

For experiment IV
1. Run exp1data.R to get simulated data sets
2. Run cv2n.R, cv3n.R to get result of CV and choose parameter for splines
3. Run exp1n.R, exp1nn.R, exp1nnn.R to get estimation results of n=400, 600, 800
4. Run exp1nIBS.R, exp1nCindex.R, exp1MAE.R to get prediction results of IBS, C-index, logrank test and MAE
5. Run result.R to analysis estimation results from exp1n.R, exp1nn.R, exp1nnn.R

For Application CGD data set

1. Run 1.cgdz2.cv.R to get CV of recurrent event model of all varying coefficient
2. Run 2.cgdz2.R to get estimation and analysis of recurrent event model of all varying coefficient
3. Run 3.cgdzpx1cv.R to get CV of recurrent event model of selected varying coefficient
4. Run 4.cgdzpx1.R to get estimation and analysis of recurrent event model of selected varying coefficient
5. Run cgdkmfunction.R to get km estimation function plot








