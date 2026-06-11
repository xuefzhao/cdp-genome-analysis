version 1.0

workflow sample_variant_summary {
    input {
        File vcf_gz
        File vcf_gz_tbi
        String sv_base_docker
    }

    call summarize_variants {
        input:
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            sv_base_docker = sv_base_docker

    }

    output {
        File sample_summary = summarize_variants.summary_tsv
    }
}

task summarize_variants {

    input {
        File vcf_gz
        File vcf_gz_tbi
        String  sv_base_docker
    }

    command <<<
        set -euo pipefail

        bcftools query -l ~{vcf_gz} > samples.txt

        echo -e "sample\tnonref_total\thetero_ref_nonref\thomo_nonref\tcompound_het\tnull_gt\tnonref_snv\tnonref_del\tnonref_ins" > sample_summary.tsv

        while read SAMPLE
        do

            bcftools query \
                -s ${SAMPLE} \
                -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' \
                ~{vcf_gz} \
            | awk -v sample="${SAMPLE}" '

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
                ref=$3
                altstr=$4
                gt=$5

                gsub(/\|/,"/",gt)

                nalt=split(altstr,alts,",")

                if(gt ~ /\./){
                    nullgt++
                    next
                }

                split(gt,g,"/")

                a1=g[1]
                a2=g[2]

                class=""

                if((a1==0 && a2>0) || (a2==0 && a1>0)){
                    het++
                    nonref++
                    class="nonref"
                }
                else if(a1>0 && a2>0 && a1==a2){
                    hom++
                    nonref++
                    class="nonref"
                }
                else if(a1>0 && a2>0 && a1!=a2){
                    comphet++
                    nonref++
                    class="nonref"
                }

                if(class=="nonref"){

                    used[a1]=1
                    used[a2]=1

                    for(idx in used){

                        if(idx==0) continue

                        alt=alts[idx]

                        if(length(ref)==length(alt)){
                            snv++
                        }
                        else if(length(ref)>length(alt)){
                            del++
                        }
                        else{
                            ins++
                        }
                    }

                    delete used
                }
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
            }' >> sample_summary.tsv

        done < samples.txt
    >>>

    output {
        File summary_tsv = "sample_summary.tsv"
    }

    runtime {
        docker: sv_base_docker
        memory: "4G"
        cpu: 1
    }
}