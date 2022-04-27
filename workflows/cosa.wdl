version 1.0

workflow run_cosa {
	input {
		String accession
		File primer_trimmed_fastq

		File reference
		File reference_index

		String container_registry
	}

	call align {
		input:
			samplename = accession,
			primer_trimmed_fastq = primer_trimmed_fastq,
			reference = reference,
			container_registry = container_registry
	}

	call deepvariant {
		input:
			samplename = accession,
			aligned_bam = align.aligned_bam,
			aligned_bam_index = align.aligned_bam_index,
			reference = reference,
			reference_index = reference_index,
			container_registry = container_registry
	}

	call VCF_consensus_deepvariant {
		input:
			samplename = accession,
			aligned_bam = align.aligned_bam,
			aligned_bam_index = align.aligned_bam_index,
			vcf = deepvariant.vcf,
			reference = reference,
			container_registry = container_registry
	}

	output {
		File aligned_bam = align.aligned_bam
		File aligned_bam_index = align.aligned_bam_index
		File high_quality_variants = VCF_consensus_deepvariant.high_quality_variants
		File high_quality_variants_index = VCF_consensus_deepvariant.high_quality_variants_index
		File consensus_fa = VCF_consensus_deepvariant.consensus_fa
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
	}
}

task align {
	input {
		String samplename
		File primer_trimmed_fastq
		File reference

		String container_registry
	}

	command {
		minimap2 \
			-a ~{reference} \
			~{primer_trimmed_fastq} \
			> ~{samplename}.aligned.sam

		samtools view \
			-bS ~{samplename}.aligned.sam \
		| samtools sort \
			> ~{samplename}.aligned.bam

		samtools index ~{samplename}.aligned.bam
	}

	output {
		File aligned_bam = "~{samplename}.aligned.bam"
		File aligned_bam_index = "~{samplename}.aligned.bam.bai"
	}

	runtime {
		docker: "~{container_registry}/cosa:c0a2fa8"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk 100 HDD"
	}
}

task deepvariant {
	input {
		String samplename
		File aligned_bam
		File aligned_bam_index
		File reference
		File reference_index

		String container_registry
	}

	Int threads = 16

	command {
		run_deepvariant \
			--model_type PACBIO \
			--ref ~{reference} \
			--reads ~{aligned_bam} \
			--output_vcf ~{samplename}.vcf \
			--num_shards ~{threads}
	}

	output {
		File vcf = "~{samplename}.vcf"
	}

	runtime {
		docker: "~{container_registry}/deepvariant:1.1.0"
		cpu: threads
		memory: "32 GB"
		disks: "local-disk 100 HDD"
	}
}

task VCF_consensus_deepvariant {
	input {
		String samplename
		File aligned_bam
		File aligned_bam_index
		File vcf
		File reference

		String container_registry
	}

	command {
		samtools mpileup \
			--min-BQ 1 \
			-f ~{reference} \
			-s ~{aligned_bam} \
			> ~{samplename}.aligned.bam.mpileup

		samtools depth \
			-q 0 \
			-Q 0 \
			~{aligned_bam} \
			> ~{samplename}.aligned.bam.depth

		VCFCons.py \
			~{reference} \
			~{samplename} \
			-c 4 \
			-f 0.5 \
			--vcf_type deepvariant \
			-q 0 \
			--input_depth ~{samplename}.aligned.bam.depth \
			--input_vcf ~{vcf}

		# change fasta header name to just be the samplename
		sed "1s/.*/>~{samplename}/" ~{samplename}.vcfcons.fasta > ~{samplename}.fasta

		minimap2 \
			-a ~{reference} \
			~{samplename}.vcfcons.frag.fasta \
			> ~{samplename}.vcfcons.frag.fasta.sam

		samtools view \
			-bS ~{samplename}.vcfcons.frag.fasta.sam \
			> ~{samplename}.vcfcons.frag.fasta.bam

		samtools sort \
			~{samplename}.vcfcons.frag.fasta.bam \
			> ~{samplename}.vcfcons.frag.fasta.sorted.bam

		samtools index \
			~{samplename}.vcfcons.frag.fasta.sorted.bam

		# fix sample name and GQ field type (to be able to import via variant transforms)
		bcftools reheader --samples <(echo "~{samplename}") ~{samplename}.vcfcons.vcf -o ~{samplename}.sample_fixed.vcf
		bcftools view -h --no-version ~{samplename}.sample_fixed.vcf > header.txt
		sed -i 's/##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Conditional genotype quality">/##FORMAT=<ID=GQ,Number=1,Type=Float,Description="Conditional genotype quality">/' header.txt
		bcftools reheader -h header.txt ~{samplename}.sample_fixed.vcf -o ~{samplename}.vcf

		bgzip ~{samplename}.vcf
		tabix ~{samplename}.vcf.gz
	}

	output {
		File high_quality_variants = "~{samplename}.vcf.gz"
		File high_quality_variants_index = "~{samplename}.vcf.gz.tbi"
		File consensus_fa = "~{samplename}.fasta"
	}

	runtime {
		docker: "~{container_registry}/cosa:c0a2fa8"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk 100 HDD"
	}
}
