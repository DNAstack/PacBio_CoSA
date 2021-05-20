version 1.0

workflow cosa {
	input {
		String samplename
		File? primer_trimmed_bam
		File? primer_trimmed_fastq
		File reference
		File reference_index
	}

	if (defined(primer_trimmed_bam)) {
		call bam_to_fastq {
			input:
				samplename = samplename,
				primer_trimmed_bam = select_first([primer_trimmed_bam])
		}
	}

	File input_fastq = if (defined(primer_trimmed_fastq)) then select_first([primer_trimmed_fastq]) else select_first([bam_to_fastq.primer_trimmed_fastq])

	call align {
		input:
			samplename = samplename,
			primer_trimmed_fastq = input_fastq,
			reference = reference
	}

	call deepvariant {
		input:
			samplename = samplename,
			aligned_bam = align.aligned_bam,
			aligned_bam_index = align.aligned_bam_index,
			reference = reference,
			reference_index = reference_index
	}

	call VCF_consensus_deepvariant {
		input:
			samplename = samplename,
			aligned_bam = align.aligned_bam,
			aligned_bam_index = align.aligned_bam_index,
			vcf = deepvariant.vcf,
			reference = reference
	}

	call assign_lineage {
		input:
			samplename = samplename,
			consensus_sequence = VCF_consensus_deepvariant.consensus_sequence
	}

	output {
		File aligned_bam = align.aligned_bam
		File aligned_bam_index = align.aligned_bam_index
		File variants = VCF_consensus_deepvariant.high_quality_variants
		File variants_index = VCF_consensus_deepvariant.high_quality_variants_index
		File consensus_sequence = VCF_consensus_deepvariant.consensus_sequence
		File consensus_sequence_fragments = VCF_consensus_deepvariant.consensus_sequence_fragments
		File lineage_metadata = assign_lineage.lineage_metadata
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
	}
}

task bam_to_fastq {
	input {
		String samplename
		File primer_trimmed_bam
	}

	Int disk_size = ceil(size(primer_trimmed_bam, "GB") * 2 + 50)

	command {
		bamtools convert \
			-format fastq \
			-in ~{primer_trimmed_bam} \
			> ~{samplename}.primer_trimmed.fastq
	}

	output {
		File primer_trimmed_fastq = "~{samplename}.primer_trimmed.fastq"
	}

	runtime {
		docker: "dnastack/cosa:c7dc46d"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " HDD"
	}
}

task align {
	input {
		String samplename
		File primer_trimmed_fastq
		File reference
	}

	Int disk_size = ceil(size(primer_trimmed_fastq, "GB") * 4 + 50)

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
		docker: "dnastack/cosa:c7dc46d"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " HDD"
	}
}

task deepvariant {
	input {
		String samplename
		File aligned_bam
		File aligned_bam_index
		File reference
		File reference_index
	}

	Int threads = 16
	Int disk_size = ceil(size(aligned_bam, "GB") * 4 + 50)

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
		docker: "google/deepvariant:1.1.0"
		cpu: threads
		memory: "32 GB"
		disks: "local-disk " + disk_size + " HDD"
	}
}

task VCF_consensus_deepvariant {
	input {
		String samplename
		File aligned_bam
		File aligned_bam_index
		File vcf
		File reference
	}

	Int disk_size = ceil(size(aligned_bam, "GB") * 4 + 50)

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

		bgzip ~{samplename}.vcfcons.vcf
		tabix ~{samplename}.vcfcons.vcf.gz
	}

	output {
		File high_quality_variants = "~{samplename}.vcfcons.vcf.gz"
		File high_quality_variants_index = "~{samplename}.vcfcons.vcf.gz.tbi"
		File consensus_sequence = "~{samplename}.vcfcons.fasta"
		File consensus_sequence_fragments = "~{samplename}.vcfcons.frag.fasta"
		File fragments_bam = "~{samplename}.vcfcons.frag.fasta.sorted.bam"
		File fragments_bam_index = "~{samplename}.vcfcons.frag.fasta.sorted.bam.bai"
	}

	runtime {
		docker: "dnastack/cosa:c7dc46d"
		cpu: 2
		memory: "7.5 GB"
		disks: "local-disk " + disk_size + " HDD"
	}
}

task assign_lineage {
	input {
		String samplename
		File consensus_sequence
	}

	Int threads = 1

	command {
		pangolin \
			--threads ~{threads} \
			--outfile ~{samplename}.lineage.csv \
			~{consensus_sequence}
	}

	output {
		File lineage_metadata = "~{samplename}.lineage.csv"
	}

	runtime {
		docker: "dnastack/pangolin:a1f8a3a"
		cpu: threads
		memory: "3.75 GB"
		disks: "local-disk 50 HDD"
	}
}
