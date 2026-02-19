nextflow.preview.output = true

// NUCLEAZE filters the reads for those containing k-mers of the target sequence
// Input: Paired-end reads (zstd compressed and interleaved), target sequence/k-mers of target sequence, k used for Nucleaze
//    First input argument is a tuple (list of sample_divs, list of paths)
// Output: Filtered forward and reverse reads that match the reference
process NUCLEAZE {
  label "medium"
  label "Nucleaze"
  input:
    tuple val(sample_divs), path(reads_files)
    path(ref_fasta_path)
    val(k)
  output:
    path("*_1.fastq"), emit: fwd_reads
    path("*_2.fastq"), emit: rev_reads
  script:
    """
    # Convert space-separated lists to arrays
    IFS=' ' read -ra samples <<< "${sample_divs.join(' ')}"
    IFS=' ' read -ra files <<< "${reads_files.join(' ')}"

    # Process each file with its corresponding sample_div
    for i in "\${!samples[@]}"; do
      sample_div="\${samples[\$i]}"
      reads_file="\${files[\$i]}"

      # Run nucleaze with equivalent parameters to former BBDuk usage:
      #   rcomp=t -> --canonical
      #   minkmerhits=1 -> --minhits 1
      #   mm=f -> (default, exact matching)
      #   interleaved=t -> --interinput
      #   ordered=t -> --order
      zstdcat "\${reads_file}" | nucleaze \\
        --in - \\
        --ref ${ref_fasta_path} \\
        --outm "\${sample_div}_1.fastq" \\
        --outm2 "\${sample_div}_2.fastq" \\
        --outu /dev/null \\
        --outu2 /dev/null \\
        --k ${k} \\
        --canonical \\
        --minhits 1 \\
        --interinput \\
        --order \\
        --threads ${task.cpus}
    done
    """
}

// CONCAT_READS concatenates several forward and reverse reads into a single forward and reverse read file
// Read headers are annotated with the basename of the file they came from
process CONCAT_READS{
  label "small"
  label "coreutils"
  input:
    path(fwd_reads)
    path(rev_reads)
  output:
    path("reads_1.fastq"), emit: final_fwd_read
    path("reads_2.fastq"), emit: final_rev_read
  script:
    """
    # Check that forward and reverse read counts match
    [ \$(echo ${fwd_reads} | wc -w) -eq \$(echo ${rev_reads} | wc -w) ] || {
      echo "Error: Unequal number of forward (\$(echo ${fwd_reads} | wc -w)) and reverse (\$(echo ${rev_reads} | wc -w)) read files"
      exit 1
    }

    # Process forward reads
    for fwd in ${fwd_reads}; do
      filename=\$(basename "\$fwd" _1.fastq)
      awk -v fname="\$filename" 'NR%4==1 {print \$0 " " fname; next} {print}' "\$fwd"
    done > "reads_1.fastq"

    # Process reverse reads
    for rev in ${rev_reads}; do
      filename=\$(basename "\$rev" _2.fastq)
      awk -v fname="\$filename" 'NR%4==1 {print \$0 " " fname; next} {print}' "\$rev"
    done > "reads_2.fastq"
    """
}

workflow {
  main:
    // Load input files and create channel
    reads = Channel.fromPath(params.s3_files)
      .splitCsv(header: false, sep: "\t")
      .map { row -> tuple(row[0], row[1]) }

    // Batch inputs
    batch_size = params.batch_size ?: 10
    batched_reads = reads
      .collate(batch_size)
      .map { batch ->
        def sample_divs = batch.collect { it[0] }
        def files = batch.collect { file(it[1]) }
        tuple(sample_divs, files)
      }

    // Process batches
    nucleaze_results = NUCLEAZE(batched_reads, params.ref_fasta_path, params.kmer)

    // Collect and filter outputs, removing empty (no-hit) files
    fwd_reads = nucleaze_results.fwd_reads
      .flatten()
      .filter { it.size() > 0 }
      .toSortedList { a, b ->
        a.name.tokenize('_')[0] <=> b.name.tokenize('_')[0]
      }
    rev_reads = nucleaze_results.rev_reads
      .flatten()
      .filter { it.size() > 0 }
      .toSortedList { a, b ->
        a.name.tokenize('_')[0] <=> b.name.tokenize('_')[0]
      }
    CONCAT_READS(fwd_reads, rev_reads)

  publish:
    fwd_reads = CONCAT_READS.out.final_fwd_read
    rev_reads = CONCAT_READS.out.final_rev_read
}

output {
  fwd_reads {
    path "results"
  }
  rev_reads {
    path "results"
  }
}
