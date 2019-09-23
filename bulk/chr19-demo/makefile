LOCAL_DATA_DIRECTORY=/Users/paul/github/ChIPseqMotifMatch/bulk/chr19-demo
BAM_FILE=ctcf_chr19.bam
NAME=demo

peaks:
	docker run -v $(LOCAL_DATA_DIRECTORY):/data/ \
               fooliu/macs2 callpeak \
              -t /data/$(BAM_FILE) \
              -n $(NAME) \
              --outdir /data/
