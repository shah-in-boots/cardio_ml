annotators folder: Contains latest code for developing ECG annotators

ludb_set.RData: .RData file containing a list of all 200 LUDB samples. Contains 12 lead of ECG signal and accompanying 12 annotation leads

create_ludb_set.R: script to create the ludb_set.RData. Largely for reference. 

bilstm_cnn_uih.R: script to train the model
	- CNN --> biLSTM --> output
	1. prepare training dataset (both LUDB and supplemental UIH data)
	2. build model
	3. train model
	4. test model
	5. save, add to model log

annotator_prep_functions.R: Contains functions for creating training and testing signal and annotation matrices, ready for ML input
    prep_ludb()	
    	Input: 
        	preferred lead (integer of value 1 thru 12)
        	annotator style (see file, may leave blank)
        	split (% samples in training set, may leave blank)
		dilation_range: percentage to dilate/compress
		max_noise: max amplitude of noise to add, as a percentage of the total ECG amplitude
		rounds: number of times to duplicate the training set 
		numer_of_derivs: number of derivatives to include (usually 2)
   	 Output:
        	list of:
            	training/testing ECG signal / annotations, in correct format for ML input
            	training/testing samples: LUDB samples used in each set. Largely for reference, usually not needed
