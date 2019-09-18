SAMPLE=GSM749704
BAMFILE=GSM749704_hg19_wgEncodeUwTfbsGm12878CtcfStdAlnRep1.bam
WORKDIR=$(HOME)/github/ChIPseqMotifMatch/bulk
default:
	@echo targets: peaks motifs merge analyze

peaks: $(BAMFILE) narrowPeaks broadPeaks

narrowPeaks:
	docker run -v $(WORKDIR)/$(SAMPLE):/data/ \
               fooliu/macs2 callpeak \
              -t /data/$(BAMFILE) \
              -n $(SAMPLE) \
              --outdir /data/

broadPeaks:
	docker run -v $(WORKDIR)/$(SAMPLE):/data/ \
               fooliu/macs2 callpeak \
              -t /data/$(BAMFILE) \
              -n $(SAMPLE) \
              --broad \
              --outdir /data/

