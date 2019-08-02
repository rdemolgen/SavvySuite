Ignored_SRC=SavvyVcfHomozygosity
SRC = CoverageBinner FishersExactTest Interval PrepareLinkageData SavvyCNV SavvyHomozygosity SavvySharedHaplotypes VariantArray Main AbstractApplication
ifndef GATK
$(error $$GATK compile using 'make GATK=/path/to/gatk/3.*/GenomeAnalysisTK.jar)
endif

savvysuite.jar : $(addprefix src/savvy/,$(addsuffix .java, $(SRC)))
	rm -rf tmp
	mkdir tmp
	echo "Manifest-Version: 1.0" > tmp.mf
	echo "Main-Class: savvy.Main" >> tmp.mf
	echo "Class-Path: $(realpath $(GATK))" | fold -w 71 | awk '{printf("%s%s\n",(NR==1?"": " "),$$0);}' >>  tmp.mf
	javac -d tmp -sourcepath src -cp "$(GATK):." $^
	jar cmvf tmp.mf $@ -C tmp .
	rm -rf tmp tmp.mf
	@echo "executable jar is: $@"


clean:
	rm -rf tmp savvysuite.jar
