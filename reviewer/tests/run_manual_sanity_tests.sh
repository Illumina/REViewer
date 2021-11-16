#!/bin/bash
TEST_OUTPUTS=sanity_test_outputs.html

for bamlet_file in inputs/bamlets/*.bam; do
  filename=$(basename -- ${bamlet_file})
  basename=${filename%.bam}
  locus=${basename%_*}
  sample=${basename#*_}

  echo "Running REViewer on ${bamlet_file}"
  reviewer \
    --reads $bamlet_file \
    --vcf inputs/vcfs/${sample}.vcf \
    --reference inputs/genomes/HG38_chr16.fa \
    --catalog inputs/catalogs/stranger_variant_catalog_hg38_chr16.json \
    --locus ${locus} \
    --output-prefix outputs/images/${basename}
done

python process_sanity_test_outputs.py outputs/images/ frozen_outputs/images/

if [[ "$OSTYPE" == "darwin"* ]]; then
  echo "Detected OS type 'dwarwin'."
  open $TEST_OUTPUTS
  if [ $? -eq 0 ]; then
    OUTPUT_DISPLAY_ERROR=0
  else
    echo "Failed to view test outputs."
    OUTPUT_DISPLAY_ERROR=1
  fi
elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
  echo "Detected OS type 'linux-gnu'"
  xdg-open $TEST_OUTPUTS
  if [ $? -eq 0 ]; then
    OUTPUT_DISPLAY_ERROR=0
  else
    echo "Failed to view test outputs."
    OUTPUT_DISPLAY_ERROR=1
  fi
else
  echo "Unknown OS type."
  OUTPUT_DISPLAY_ERROR=1
fi

if [ OUTPUT_DISPLAY_ERROR ]; then
  echo "To verify sanity tests, please view ${TEST_OUTPUTS} in a browser."
fi