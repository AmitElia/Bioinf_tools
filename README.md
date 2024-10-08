Variantcalling Makefile:
  make init run ACC=${ACC} SRR=${SRR} N=${N}"
	
	make init runbatch ACC=${ACC} N=${N} PRJN=${PRJN} NR=${NR}"
	
	make annotate ACC=${ACC} PRJN=${PRJN}"
	
	make visualize ACC=${ACC} SRR=${SRR}"
	
	init - initialize the genome."
	
	run - get fastq, trim, align, and call variants."
	
	runbatch - fastq, trim, align, and call variants based on PRJN bioproject number."
	
	annotate - creates custom annotation based on ACC number."
	
	visualize - visualize alignment using IGV"
	
	clean - remove all files created."
	
