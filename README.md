# Analysing Cancer Data

This project aims to explore and analyse cancer data for educational purposes.

## Download data

Install AWS CLI from: https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html

Run the following command from the root directory of this project:

```bash
aws s3 cp s3://icgc25k-open/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz data --endpoint-url https://object.genomeinformatics.org --no-sign-request
```