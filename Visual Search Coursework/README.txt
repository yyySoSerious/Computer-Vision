Source code is found in the cwork_basecode_2012 folder
-----------------------------------------------------------------------------------------------------
=====> To save a descriptor (globalRGB, spatial Grid , eoh, eoh with colour):
	--->  Go to cvpr_computedescriptors.m
	--->  From line 25, uncomment the descriptor you want to save and press play

=====> To save the visual words descriptor:
	---> First, Save the sift descriptors by going to computeSiftDescriptors and pressing play.
	---> Then, go to visualWords.m and pressing play (although, I will not recommend it because 	     it is really slow at the current configuration). If you just want to see it running 
	     without caring about the performance of the system, you can comment out the line 2 in 
 	     siftDescriptors.m, which will increase the speed of the algorithm, although it is still 
  	     not very fast because it is processing over 170000 descriptors. You can also increase 	     speed by reducing the max number of iterations of the k-means algorithm by changing the 	     last argument on line 35. It is currently set to 1000.

=====> To run a visual search on a system, plot pr-curve and visualize result (The descriptor must be 
       saved beforehand):
	---> Go to cvpr_visualsearch.m
	---> From line 30, uncomment the descriptor, you want to run
	---> If you want to change the distance measure (excluding Mahalanobis), go to line 57, and 	     uncomment the plotPRCurve function that matches your selected choice.
	---> if you want to change to Mahalanobis distance, go to ComputeDescriptorsWithPCA , choose 	     a descriptor from line 7, by uncommenting the chosen line. then press play. When done, 	     go back to cvpr_visualsearch.m, uncomment your selected descriptor from line 30, then 	     comment out all the other plotPRCurve functions from line 57 and uncomment only lines 	     64, 65 and 66. Finally, press play.

=====> To use Mahalanobis as a cost metric:
	---> Go to computeDescriptorsWithPCA
	---> Choose the descriptor type from line 7, then press play.
	---> go to cvpr_visualsearch.m
	---> Select you preferred descriptor on line 30.
	---> comment out the plotPRCurve functions from line 57 and uncomment lines 71, 72 and 73, 	     then press play.

