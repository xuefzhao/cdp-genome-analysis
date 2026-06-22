version 1.0
import "../utils/Structs.wdl" as Structs

workflow vcf_to_bed_and_variant_class_counts {
  input {
    Array[File] vcfs
    File vcf_to_bed_script
    File count_variant_classes_script
    File merge_counts_script
    File plot_hist_script
    File sample_tsv
    File batch_color_tsv
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_vcf_to_bed
    RuntimeAttr? runtime_attr_count_variant_classes_per_contig
    RuntimeAttr? runtime_attr_concat_beds
    RuntimeAttr? runtime_attr_merge_synonymous_counts
    RuntimeAttr? runtime_attr_merge_lof_counts
    RuntimeAttr? runtime_attr_merge_missense_counts
    RuntimeAttr? runtime_attr_merge_others_counts
    RuntimeAttr? runtime_attr_plot_variant_histograms
  }

  scatter (vcf in vcfs) {
    call vcf_to_bed {
      input:
        vcf = vcf,
        script = vcf_to_bed_script,
        docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_vcf_to_bed
    }

    call count_variant_classes_per_contig {
      input:
        bed_tsv = vcf_to_bed.bed_tsv,
        script = count_variant_classes_script,
        docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_count_variant_classes_per_contig
    }
  }

  call concat_beds {
    input:
      bed_files = vcf_to_bed.bed_tsv,
      docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_concat_beds
  }

  call merge_synonymous_counts {
    input:
      count_tables = count_variant_classes_per_contig.synonymous_count_table,
      script = merge_counts_script,
      docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_merge_synonymous_counts
  }

  call merge_lof_counts {
    input:
      count_tables = count_variant_classes_per_contig.lof_count_table,
      script = merge_counts_script,
      docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_merge_lof_counts
  }

  call merge_missense_counts {
    input:
      count_tables = count_variant_classes_per_contig.missense_count_table,
      script = merge_counts_script,
      docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_merge_missense_counts
  }

  call merge_others_counts {
    input:
      count_tables = count_variant_classes_per_contig.others_count_table,
      script = merge_counts_script,
      docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_merge_others_counts
  }

  call plot_variant_histograms {
    input:
      synonymous_table = merge_synonymous_counts.merged_count_table,
      lof_table = merge_lof_counts.merged_count_table,
      missense_table = merge_missense_counts.merged_count_table,
      others_table = merge_others_counts.merged_count_table,
      sample_tsv = sample_tsv,
      batch_color_tsv = batch_color_tsv,
      script = plot_hist_script,
      docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_plot_variant_histograms
  }

  output {
    File united_bed_gz = concat_beds.united_bed_gz
    File united_bed_gz_tbi = concat_beds.united_bed_gz_tbi

    File synonymous_counts_per_sample = merge_synonymous_counts.merged_count_table
    File lof_counts_per_sample = merge_lof_counts.merged_count_table
    File missense_counts_per_sample = merge_missense_counts.merged_count_table
    File others_counts_per_sample = merge_others_counts.merged_count_table

    File synonymous_histogram_png = plot_variant_histograms.synonymous_histogram_png
    File lof_histogram_png = plot_variant_histograms.lof_histogram_png
    File missense_histogram_png = plot_variant_histograms.missense_histogram_png
    File others_histogram_png = plot_variant_histograms.others_histogram_png
  }
}

task vcf_to_bed {
  input {
    File vcf
    File script
    String docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: 5 * ceil(size(vcf, "GB")) + 10,
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    python3 "~{script}" --vcf "~{vcf}" --out out.bed.tsv
  >>>

  output {
    File bed_tsv = "out.bed.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: docker
  }
}

task count_variant_classes_per_contig {
  input {
    File bed_tsv
    File script
    String docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 5 * ceil(size(bed_tsv, "GB")) + 10,
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    python3 "~{script}" --bed "~{bed_tsv}" --out-prefix class_counts
  >>>

  output {
    File synonymous_count_table = "class_counts.synonymous.stat.tsv"
    File lof_count_table = "class_counts.lof.stat.tsv"
    File missense_count_table = "class_counts.missense.stat.tsv"
    File others_count_table = "class_counts.others.stat.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: docker
  }
}

task concat_beds {
  input {
    Array[File] bed_files
    String docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: 5 * ceil(size(bed_files, "GB")) + 10,
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    if [ ~{length(bed_files)} -eq 0 ]; then
      echo "No BED files provided." >&2
      exit 1
    fi

    head -n 1 "~{bed_files[0]}" > united.annotated.bed
    for f in ~{sep=' ' bed_files}; do
      tail -n +2 "$f"
    done | sort -t$'\t' -k1,1 -k2,2n >> united.annotated.bed

    bgzip -c united.annotated.bed > united.annotated.bed.gz
    tabix -S 1 -p bed united.annotated.bed.gz
  >>>

  output {
    File united_bed = "united.annotated.bed"
    File united_bed_gz = "united.annotated.bed.gz"
    File united_bed_gz_tbi = "united.annotated.bed.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: docker
  }
}

task merge_synonymous_counts {
  input {
    Array[File] count_tables
    File script
    String docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 5 * ceil(size(count_tables, "GB")) + 10,
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    python3 "~{script}" --tables ~{sep=' ' count_tables} --out synonymous_counts_per_sample.tsv
  >>>
  output {
    File merged_count_table = "synonymous_counts_per_sample.tsv"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: docker
  }
}

task merge_lof_counts {
  input {
    Array[File] count_tables
    File script
    String docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 5 * ceil(size(count_tables, "GB")) + 10,
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    python3 "~{script}" --tables ~{sep=' ' count_tables} --out lof_counts_per_sample.tsv
  >>>
  output {
    File merged_count_table = "lof_counts_per_sample.tsv"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: docker
  }
}

task merge_missense_counts {
  input {
    Array[File] count_tables
    File script
    String docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 5 * ceil(size(count_tables, "GB")) + 10,
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    python3 "~{script}" --tables ~{sep=' ' count_tables} --out missense_counts_per_sample.tsv
  >>>
  output {
    File merged_count_table = "missense_counts_per_sample.tsv"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: docker
  }
}

task merge_others_counts {
  input {
    Array[File] count_tables
    File script
    String docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 5 * ceil(size(count_tables, "GB")) + 10,
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    python3 "~{script}" --tables ~{sep=' ' count_tables} --out others_counts_per_sample.tsv
  >>>
  output {
    File merged_count_table = "others_counts_per_sample.tsv"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: docker
  }
}

task plot_variant_histograms {
  input {
    File synonymous_table
    File lof_table
    File missense_table
    File others_table
    File sample_tsv
    File batch_color_tsv
    File script
    String docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: 5 * ceil(size([synonymous_table, lof_table, missense_table, others_table, sample_tsv, batch_color_tsv], "GB")) + 10,
    boot_disk_gb: 10,
    preemptible_tries: 2,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    Rscript "~{script}" \
      --synonymous "~{synonymous_table}" \
      --lof "~{lof_table}" \
      --missense "~{missense_table}" \
      --others "~{others_table}" \
      --sample_tsv "~{sample_tsv}" \
      --batch_color_tsv "~{batch_color_tsv}" \
      --out_prefix "variant_count_per_sample"
  >>>
  output {
    File synonymous_histogram_png = "variant_count_per_sample.synonymous.hist.png"
    File lof_histogram_png = "variant_count_per_sample.lof.hist.png"
    File missense_histogram_png = "variant_count_per_sample.missense.hist.png"
    File others_histogram_png = "variant_count_per_sample.others.hist.png"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: docker
  }
}
