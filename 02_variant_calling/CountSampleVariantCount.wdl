version 1.0

workflow CountSampleVariantCount {

    input {
        File vcf_gz
        File vcf_gz_tbi
        String sv_base_pipeline
    }

    call extract_samples {
        input:
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            sv_base_pipeline = sv_base_pipeline
    }

    scatter(sample in extract_samples.samples) {

        call summarize_sample {
            input:
                vcf_gz = vcf_gz,
                vcf_gz_tbi = vcf_gz_tbi,
                sv_base_pipeline = sv_base_pipeline,
                sample = sample
        }
    }

    call merge_results {
        input:
            sample_tables = summarize_sample.sample_summary,
            sv_base_pipeline = sv_base_pipeline
    }

    output {
        File sample_summary = merge_results.summary_tsv
    }
}


task extract_samples {

    input {
        File vcf_gz
        File vcf_gz_tbi
        String sv_base_pipeline
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf_gz} > samples.txt
    >>>

    output {
        Array[String] samples = read_lines("samples.txt")
    }

    runtime {
        docker: sv_base_pipeline
        memory: "1G"
    }
}

task summarize_sample {

    input {
        String sample
        File vcf_gz
        File vcf_gz_tbi
        String sv_base_pipeline
    }

    command <<<
        set -euo pipefail

        echo -e "sample\tnonref_total\thetero_ref_nonref\thomo_nonref\tcompound_het\tnull_gt\tnonref_snv\tnonref_del\tnonref_ins" > result.tsv

        bcftools query \
            -s ~{sample} \
            -f '%REF\t%ALT[\t%GT]\n' \
            ~{vcf_gz} \
        | awk -v sample="~{sample}" '

        BEGIN{
            nonref=0
            het=0
            hom=0
            comphet=0
            nullgt=0
            snv=0
            del=0
            ins=0
        }

        {
            ref=$1
            altstr=$2
            gt=$3

            gsub(/\|/,"/",gt)

            if(gt ~ /\./){
                nullgt++
                next
            }

            split(gt,g,"/")
            a1=g[1]
            a2=g[2]

            split(altstr,alts,",")

            if((a1==0 && a2>0) || (a2==0 && a1>0)){
                het++
                nonref++
            }
            else if(a1>0 && a2>0 && a1==a2){
                hom++
                nonref++
            }
            else if(a1>0 && a2>0 && a1!=a2){
                comphet++
                nonref++
            }
            else{
                next
            }

            used[a1]=1
            used[a2]=1

            for(idx in used){

                if(idx==0) continue

                alt=alts[idx]

                if(length(ref)==length(alt))
                    snv++
                else if(length(ref)>length(alt))
                    del++
                else
                    ins++
            }

            delete used
        }

        END{
            print sample "\t" \
                  nonref "\t" \
                  het "\t" \
                  hom "\t" \
                  comphet "\t" \
                  nullgt "\t" \
                  snv "\t" \
                  del "\t" \
                  ins
        }' >> result.tsv
    >>>

    output {
        File sample_summary = "result.tsv"
    }

    runtime {
        docker: sv_base_pipeline
        memory: "2G"
    }
}

task merge_results {

    input {
        Array[File] sample_tables
        String sv_base_pipeline
    }

    command <<<
        set -euo pipefail

        head -n 1 ~{sample_tables[0]} > summary.tsv

        for f in ~{sep=' ' sample_tables}
        do
            tail -n +2 $f >> summary.tsv
        done
    >>>

    output {
        File summary_tsv = "summary.tsv"
    }

    runtime {
        docker: sv_base_pipeline
        memory: "1G"
    }
}

